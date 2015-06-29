#include <iostream>
#include <chrono>
#include <cassert>
#include <functional>
#include <cstring>
#include "parse_arguments.h"
#include "plummer.h"

template<typename T>
T square(T x) { return x * x ; }

template<typename T>
T rsqrt(T x) { return T{1}/sqrt(x) ; }


template<typename Real>
struct Hermite4T
{
  using real_type = Real;
  static constexpr auto PP_FLOPS = 44;

  const int nbody_;
  Real eta_, eps_;
  
  Real eps2_;
  Real dt_;
  Real time_;


  struct Particle
  {
    Real mass;
    Real posx, posy, posz;
    Real velx, vely, velz;
  };

  struct Force
  {
    Real accx{0}, accy{0}, accz{0};
    Real jrkx{0}, jrky{0}, jrkz{0};
    Real gpot{0};
  };

  struct Predictor
  {
    Real mass;
    Real posx, posy, posz;
    Real velx, vely, velz;
    Predictor() {} 
    Predictor(const Real dt, const Particle& p, const Force &f) noexcept
    {
      const auto dt2 = Real{1.0/2.0}*dt;
      const auto dt3 = Real{1.0/3.0}*dt;
      mass = p.mass;
      posx = p.posx + dt*(p.velx + dt2*(f.accx + dt3*f.jrkx));
      posy = p.posy + dt*(p.vely + dt2*(f.accy + dt3*f.jrky));
      posz = p.posz + dt*(p.velz + dt2*(f.accz + dt3*f.jrkz));
      velx = p.velx + dt*(f.accx + dt2* f.jrkx);
      vely = p.vely + dt*(f.accy + dt2* f.jrky);
      velz = p.velz + dt*(f.accz + dt2* f.jrkz);
    }
  };

  void set_eps(const Real eps) { eps_ = eps; }
  void set_eta(const Real eta) { eta_ = eta; }

  std::vector<Particle> ptcl_vec_;
  std::vector<Predictor> pred_vec_;
  std::vector<Force> force_vec_;
  
  Real time() const { return time_;}
  Real dt() const {return dt_;}

  Hermite4T(const int nbody, const Real Q = 1.0, const Real eta = 0.1) : 
    nbody_(nbody), eta_(eta), eps_(4.0/nbody)
  {
    eps2_ = eps_*eps_;
    time_ = 0;

    ptcl_vec_.resize(nbody_);
    pred_vec_.resize(nbody_);
    force_vec_.resize(nbody_);

    using namespace std;

    Plummer model(nbody_);

    // Scale velocity by "Q"-virial factor
    // Q=1 : maintain virial equlibirum
    // Q<1 : cold system (constracts)
    // Q>1 : hot system  (expands)
#ifdef _OPENMP
#pragma omp parallel for schedule(rtune)
#endif
    for (int i = 0; i < nbody_; ++i)
    {
      ptcl_vec_[i] = Particle{
          (Real)model.mass[i], 
          (Real)model.pos[i].x,   (Real)model.pos[i].y,   (Real)model.pos[i].z,
        Q*(Real)model.vel[i].x, Q*(Real)model.vel[i].y, Q*(Real)model.vel[i].z
      };
    }

    dt_ = Real{1.0e-4};
  }

  Force pp_force(const Force& f, const Predictor &pi, const Predictor &pj)
  {
    const auto dx = pj.posx - pi.posx;
    const auto dy = pj.posy - pi.posy;
    const auto dz = pj.posz - pi.posz;
    const auto ds2 = dx*dx + dy*dy + dz*dz + eps2_;

    const auto  inv_ds  = rsqrt(ds2);
    const auto minv_ds  = pj.mass;
    const auto  inv_ds2 = inv_ds * inv_ds;
    const auto minv_ds3 = inv_ds2*minv_ds;

    Force fi;

    fi.accx = f.accx + minv_ds3*dx;
    fi.accy = f.accy + minv_ds3*dy;
    fi.accz = f.accz + minv_ds3*dz;
    fi.gpot = f.gpot - minv_ds;

    const auto dvx = pj.velx - pi.velx;
    const auto dvy = pj.vely - pi.vely;
    const auto dvz = pj.velz - pi.velz;

    const auto rv = dx*dvx + dy*dvy + dz*dvz;

    const auto Jij = Real{-3}*(rv*inv_ds2*minv_ds3);

    fi.jrkx = f.jrkx + minv_ds3*dvx + Jij*dx;
    fi.jrky = f.jrky + minv_ds3*dvy + Jij*dy;
    fi.jrkz = f.jrkz + minv_ds3*dvz + Jij*dz;

    return fi;
  }

  std::vector<Force> compute_forces(bool no_predictor = false)
  {
    if (no_predictor)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
      for (int i = 0; i < nbody_; i++)
        pred_vec_[i] = Predictor(0, ptcl_vec_[i], force_vec_[i]);
    }

    std::vector<Force> forces(nbody_);
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif

    for (int i = 0; i < nbody_; i++)
    {
      const auto& pi = pred_vec_[i];
      for (int j = 0; j < nbody_; j++)
      {
        forces[i] = pp_force(forces[i], pi, pred_vec_[j]);
      }
    }
    return forces;
  }

  void step()
  {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
    for (int i = 0; i < nbody_; ++i)
      pred_vec_[i] = Predictor(dt_, ptcl_vec_[i], force_vec_[i]);

    const auto&& forces_new = compute_forces();

    const auto dt   = dt_;
    const auto h    = Real{0.5}*dt;
    const auto hinv = Real{1}/h;
    const auto f1   = Real{0.5}*hinv*hinv;
    const auto f2   = Real{3.0}*hinv*f1;

    const auto dt2 = dt *dt * Real{1.0/2.0};
    const auto dt3 = dt2*dt * Real{1.0/3.0};
    const auto dt4 = dt3*dt * Real{1.0/4.0};
    const auto dt5 = dt4*dt * Real{1.0/4.0};

    Real dt_min = HUGE;
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) reduction(min:dt_min)
#endif
    for (int i = 0; i < nbody_; ++i)
    {
      /* interpolate snap & crackle */
      const auto Amx = forces_new[i].accx - force_vec_[i].accx;
      const auto Amy = forces_new[i].accy - force_vec_[i].accy;
      const auto Amz = forces_new[i].accz - force_vec_[i].accz;

      const auto Jmx = h*(forces_new[i].jrkx - force_vec_[i].jrkx);
      const auto Jmy = h*(forces_new[i].jrky - force_vec_[i].jrky);
      const auto Jmz = h*(forces_new[i].jrkz - force_vec_[i].jrkz);
      
      const auto Jpx = h*(forces_new[i].jrkx + force_vec_[i].jrkx);
      const auto Jpy = h*(forces_new[i].jrky + force_vec_[i].jrky);
      const auto Jpz = h*(forces_new[i].jrkz + force_vec_[i].jrkz);

      auto snpx = f1 * Jmx;
      auto snpy = f1 * Jmy;
      auto snpz = f1 * Jmz;

      auto crkx = f2 * (Jpx - Amx);
      auto crky = f2 * (Jpy - Amy);
      auto crkz = f2 * (Jpz - Amz);

      snpx -= h*crkx;
      snpy -= h*crky;
      snpz -= h*crkz;

      /* correct position & velocity */

      ptcl_vec_[i].posx += dt4*snpx + dt5*crkx;
      ptcl_vec_[i].posy += dt4*snpy + dt5*crky;
      ptcl_vec_[i].posz += dt4*snpz + dt5*crkz;

      ptcl_vec_[i].velx += dt3*snpx + dt4*crkx;
      ptcl_vec_[i].vely += dt3*snpy + dt4*crky;
      ptcl_vec_[i].velz += dt3*snpz + dt4*crkz;

      /* compute new time-step */

      force_vec_[i] = forces_new[i];

      const auto s0 = square(force_vec_[i].accx) + square(force_vec_[i].accy) + square(force_vec_[i].accz);
      const auto s1 = square(force_vec_[i].jrkx) + square(force_vec_[i].jrky) + square(force_vec_[i].jrkz);
      const auto s2 = snpx*snpx + snpy*snpy + snpz*snpz;
      const auto s3 = crkx*crkx + crky*crky + crkz*crkz;

      using std::sqrt;
      using std::min;
      const auto u = sqrt(s0*s2) + s1;
      const auto l = sqrt(s1*s3) + s2;
      assert(l > 0.0);
      const auto dt_loc = eta_ * static_cast<Real>(sqrt(u/l));
      dt_min = min(dt_min, dt_loc);
    }

    time_ += dt_;
    dt_ = dt_min;
  }

  double etot()
  {
    double epot = 0, ekin = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ekin,epot)
#endif
    for (int i = 0; i < nbody_; i++)
    {
      const auto& p = ptcl_vec_[i];
      ekin += p.mass * (square(p.velx) + square(p.vely) + square(p.velz)) * Real{0.5};
      epot += Real{0.5} * p.mass * force_vec_[i].gpot;
    }

    return epot + ekin;
  }

};

template<typename Real>
void run(const int nbody, const int nstep, const Real eta, const Real Q)
{
  using namespace std;
  using Hermite4 = Hermite4T<Real>;
  Hermite4 h4(nbody,Q,eta);

  h4.step();

  using _time = chrono::system_clock;
  const auto etot0 = h4.etot();

  const auto start = _time::now();
  for (int step = 0; step < nstep; step++)
  {
    cerr << "  step= " << step << " out of " << nstep << ": t= " << h4.time() << "  dt= " << h4.dt() <<  endl;
    h4.step();
  }
  const auto end = _time::now();
  const auto elapsed = chrono::duration_cast<chrono::duration<double>>(end-start).count();

  const auto etot1 = h4.etot();

  cout << endl;
  cout << "Statistics: " << endl;
  cout << "\tEtot(0)= " << etot0 << endl;
  cout << "\tEtot(1)= " << etot1 << endl;
  cout << "\tdEtot=   " << abs(etot1 - etot0)/abs(etot0) << endl;
  cout << "\tElapsed time= " << elapsed << " seconds.\n";
  cout << "\tPerformance= " << 1.0*nstep*nbody*nbody*Hermite4::PP_FLOPS/elapsed/1e9 << "GFLOP/s\n";
}


int main(int argc, char *argv[])
{
  using namespace std;
  using namespace parse_arguments;

  auto nbodies = 1024;
  auto nsteps  = 16;
  auto eta     = 0.3;
  auto Q       = 1.0;
  auto fp32    = false;
//  auto fileName = std::string{};

  auto param_pack = pack(argc, argv,
      param("number of bodies",     nbodies, "n", "nbodies"),
      param("number of steps",      nsteps,  "s", "nsteps"),
      param("accuracy parameter",   eta,     "e", "eta"),
      param("virial ratio",         Q,       "Q", ""),
      param("toggle single precision", fp32,   "",  "fp32")
//      param("input filename" ,     fileName, 'f', "infile")
      );

  cerr << param_pack.parse_all([](std::string s) { cerr << " --- crap --- \n" << s; exit(2); });
#if 0
  if (!param_pack.parse())
  {
    cerr << param_pack.usage();
    exit(-1);
  }
  cerr << param_pack.params();
#endif

  assert(nbodies > 1);
  assert(nsteps > 0);
  assert(eta > 0);
  assert(Q > 0);

  if (fp32)
    run<float>(nbodies, nsteps, eta, Q);
  else
    run<double>(nbodies, nsteps, eta, Q);
}
