#include <iostream>
#include <chrono>
#include <cassert>
#include <functional>
#include <cstring>
#include "parse_arguments.h"
#include "plummer.h"
#include "safe_printf.h"

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

  void set_eps(const Real eps) { eps_ = eps; }
  void set_eta(const Real eta) { eta_ = eta; }

  std::vector<Particle> ptcl_vec_;
  std::vector<Force> force_vec_, force_new_;
  
  Real time() const { return time_;}
  Real dt() const {return dt_;}
  
  void pp_force(Force& f, const Particle &pi, const Particle &pj)
  {
    const auto dx = pj.posx - pi.posx;
    const auto dy = pj.posy - pi.posy;
    const auto dz = pj.posz - pi.posz;
    const auto r2 = dx*dx + dy*dy + dz*dz;
    const auto ds2 = r2 + eps2_;

    const auto  inv_ds  = rsqrt(ds2);
    const auto  inv_ds2 = inv_ds * inv_ds;
    const auto minv_ds  = inv_ds * pj.mass;
    const auto minv_ds3 = inv_ds2 * minv_ds;

    Force fi;

    f.accx += minv_ds3*dx;
    f.accy += minv_ds3*dy;
    f.accz += minv_ds3*dz;
    if (r2 > 0)
      f.gpot -= minv_ds;

    const auto dvx = pj.velx - pi.velx;
    const auto dvy = pj.vely - pi.vely;
    const auto dvz = pj.velz - pi.velz;

    const auto rv = dx*dvx + dy*dvy + dz*dvz;

    const auto Jij = Real{-3}*(rv*inv_ds2*minv_ds3);

    f.jrkx += minv_ds3*dvx + Jij*dx;
    f.jrky += minv_ds3*dvy + Jij*dy;
    f.jrkz += minv_ds3*dvz + Jij*dz;
  }

  void compute_forces(std::vector<Force> &forces)
  {
    for (int i = 0; i < nbody_; i++)
    {
      const auto& pi = ptcl_vec_[i];
      Force fi;
      for (int j = 0; j < nbody_; j++)
        pp_force(fi, pi, ptcl_vec_[j]);
      forces[i] = fi;
    }
  }

  Hermite4T(const int nbody, const Real Q = 1.0, const Real eta = 0.1) : 
    nbody_(nbody), eta_(eta), eps_(4.0/nbody)
  {
    eps2_ = eps_*eps_;
    time_ = 0;

    ptcl_vec_.resize(nbody_);
    force_vec_.resize(nbody_);
    force_new_.resize(nbody_);

    using namespace std;

    Plummer model(nbody_);

    // Scale velocity by "Q"-virial factor
    // Q=1 : maintain virial equlibirum
    // Q<1 : cold system (constracts)
    // Q>1 : hot system  (expands)
    for (int i = 0; i < nbody_; ++i)
    {
      ptcl_vec_[i] = Particle{
          (Real)model.mass[i], 
          (Real)model.pos[i].x,   (Real)model.pos[i].y,   (Real)model.pos[i].z,
        Q*(Real)model.vel[i].x, Q*(Real)model.vel[i].y, Q*(Real)model.vel[i].z
      };
    }

    dt_ = Real{1.0e-4};

    compute_forces(force_vec_);
  }


  void step()
  {
    {
      const auto dt  = dt_;
      const auto dt2 = dt*Real{1.0/2.0};
      const auto dt3 = dt*Real{1.0/3.0};

      for (int i = 0; i < nbody_; ++i)
      {
        auto& p = ptcl_vec_[i];
        const auto& f = force_vec_[i];

        p.posx += dt*(p.velx + dt2*(f.accx + dt3*f.jrkx));
        p.posy += dt*(p.vely + dt2*(f.accy + dt3*f.jrky));
        p.posz += dt*(p.velz + dt2*(f.accz + dt3*f.jrkz));

        p.velx += dt*(f.accx + dt2*f.jrkx);
        p.vely += dt*(f.accy + dt2*f.jrky);
        p.velz += dt*(f.accz + dt2*f.jrkz);
      }
    }

    compute_forces(force_new_);

    {
      const auto dt   = dt_;
      const auto h    = Real{0.5}*dt;
      const auto hinv = Real{1}/h;
      const auto f1   = Real{0.5}*hinv*hinv;
      const auto f2   = Real{3.0}*hinv*f1;

      const auto dt2 = dt *dt * Real{1.0/2.0};
      const auto dt3 = dt2*dt * Real{1.0/3.0};
      const auto dt4 = dt3*dt * Real{1.0/4.0};
      const auto dt5 = dt4*dt * Real{1.0/5.0};

      Real dt_min = HUGE;
      for (int i = 0; i < nbody_; ++i)
      {
        /* interpolate snap & crackle */
        const auto Amx = force_new_[i].accx - force_vec_[i].accx;
        const auto Amy = force_new_[i].accy - force_vec_[i].accy;
        const auto Amz = force_new_[i].accz - force_vec_[i].accz;

        const auto Jmx = h*(force_new_[i].jrkx - force_vec_[i].jrkx);
        const auto Jmy = h*(force_new_[i].jrky - force_vec_[i].jrky);
        const auto Jmz = h*(force_new_[i].jrkz - force_vec_[i].jrkz);

        const auto Jpx = h*(force_new_[i].jrkx + force_vec_[i].jrkx);
        const auto Jpy = h*(force_new_[i].jrky + force_vec_[i].jrky);
        const auto Jpz = h*(force_new_[i].jrkz + force_vec_[i].jrkz);

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

        force_vec_[i] = force_new_[i];

        const auto s0 = square(force_vec_[i].accx) + square(force_vec_[i].accy) + square(force_vec_[i].accz);
        const auto s1 = square(force_vec_[i].jrkx) + square(force_vec_[i].jrky) + square(force_vec_[i].jrkz);
        const auto s2 = snpx*snpx + snpy*snpy + snpz*snpz;
        const auto s3 = crkx*crkx + crky*crky + crkz*crkz;

        using std::sqrt;
        using std::min;
        const auto u = sqrt(s0*s2) + s1;
        const auto l = sqrt(s1*s3) + s2;
        assert(l > 0.0);
        const auto dt_loc = eta_ * sqrt(u/l);
        dt_min = min(dt_min, dt_loc);
      }

      time_ += dt_;
      dt_ = 1.0;
      while (dt_ > dt_min) dt_ *= 0.5;
    }
  }

  double etot()
  {
    double epot = 0, ekin = 0;
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
void run(const int nbody, const Real t_end, const int dstep, const Real eta, const Real Q)
{
  using namespace std;
  using Hermite4 = Hermite4T<Real>;
  Hermite4 h4(nbody,Q,eta);

  h4.step();

  auto etot0 = h4.etot();
  using _time = chrono::system_clock;

  const auto start = _time::now();
  auto step = 0;
  auto etot_pre = etot0;
  while (h4.time() < t_end)
  {
    if (step%dstep == 0)
    {
      const auto etot_post = h4.etot();
      const auto de = etot0 - etot_post;
      const auto dde = etot_post - etot_pre;
      etot_pre = etot_post;
      cerr << "  step= " << step << ": t= " << h4.time() << "  dt= " << h4.dt();
      cerr << " | de= " << de << " d(de)= " << dde << endl;
    }
    h4.step();
    step++;
  }
  const auto end = _time::now();
  const auto elapsed = chrono::duration_cast<chrono::duration<double>>(end-start).count();

  auto nstep = step;

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
  auto t_end   = 1.0/16;
  auto dstep   = 10;
  auto eta     = 0.4;
  auto Q       = 1.0;
  auto fp32    = false;

  auto param_pack = pack(argc, argv,
      param("number of bodies",     nbodies, "n", "nbodies"),
      param("integration time",     t_end,  "t", "tend"),
      param("print output every # steps" , dstep, "s", "dstep"),
      param("accuracy parameter",   eta,     "e", "eta"),
      param("virial ratio",         Q,       "Q", ""),
      param("toggle single precision", fp32,   "",  "fp32")
      );

  cerr << param_pack.parse_all();

  assert(nbodies > 1);
  assert(t_end > 0);
  assert(eta > 0);
  assert(Q > 0);


  if (fp32)
    run<float>(nbodies, t_end, dstep, eta, Q);
  else
    run<double>(nbodies, t_end, dstep, eta, Q);
}
