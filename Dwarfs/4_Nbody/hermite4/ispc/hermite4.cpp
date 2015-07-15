/*
  Copyright (c) 2014, Evghenii Gaburov
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.


   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Hermite4 N-body integrator */
/* Makino and Aarseth, 1992 */
/* http://adsabs.harvard.edu/abs/1992PASJ...44..141M and references there in*/

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <cassert>
#include <chrono>

#include "parse_arguments.h"
#include "typeReal.h"
#include "plummer.h"
#include "hermite4_ispc.h"

template<typename T>
T square(T x) { return x * x ; }

struct Hermite4
{
  static constexpr auto PP_FLOPS = 44;

  const int nbody_;
  Real eta_, eps_;
  Real eps2_;
  Real dt_;
  Real time_;

  // Particle
  std::vector<Real> mass_;
  std::vector<Real> posx_, posy_, posz_;
  std::vector<Real> velx_, vely_, velz_;

  // Force
  std::vector<Real> accx_, accy_, accz_;
  std::vector<Real> jrkx_, jrky_, jrkz_;
  std::vector<Real> gpot_;

  std::vector<Real> accx1_, accy1_, accz1_;
  std::vector<Real> jrkx1_, jrky1_, jrkz1_;
  
  void set_eps(const Real eps) { eps_ = eps; }
  void set_eta(const Real eta) { eta_ = eta; }
  Real time() const { return time_;}
  Real dt() const {return dt_;}

  Hermite4(const int nbody, const Real Q = 1.0, const Real eta = 0.1) : 
    nbody_(nbody), eta_(eta), eps_(4.0/nbody)
  {
    using namespace std;

    eps2_ = eps_*eps_;
    time_ = 0;

    mass_.resize(nbody_);
    posx_.resize(nbody_);
    posy_.resize(nbody_);
    posz_.resize(nbody_);
    velx_.resize(nbody_);
    vely_.resize(nbody_);
    velz_.resize(nbody_);

    accx_.resize(nbody_);
    accy_.resize(nbody_);
    accz_.resize(nbody_);
    jrkx_.resize(nbody_);
    jrky_.resize(nbody_);
    jrkz_.resize(nbody_);
    gpot_.resize(nbody_);

    accx1_.resize(nbody_);
    accy1_.resize(nbody_);
    accz1_.resize(nbody_);
    jrkx1_.resize(nbody_);
    jrky1_.resize(nbody_);
    jrkz1_.resize(nbody_);

    Plummer model(nbody_);
    // Scale velocity by "Q"-virial factor
    // Q=1 : maintain virial equlibirum
    // Q<1 : cold system (constracts)
    // Q>1 : hot system  (expands)
    for (int i = 0; i < nbody_; ++i)
    {
      mass_[i] = model.mass[i];
      posx_[i] = model.pos[i].x;
      posy_[i] = model.pos[i].y;
      posz_[i] = model.pos[i].z;
      velx_[i] = model.vel[i].x*Q;
      vely_[i] = model.vel[i].y*Q;
      velz_[i] = model.vel[i].z*Q;
    };

    dt_ = Real{1.0e-4};

    ispc::compute_forces(
        nbody_,
        mass_.data(),
        posx_.data(),
        posy_.data(),
        posz_.data(),
        velx_.data(),
        vely_.data(),
        velz_.data(),
        accx_.data(),
        accy_.data(),
        accz_.data(),
        jrkx_.data(),
        jrky_.data(),
        jrkz_.data(),
        gpot_.data(),
        eps2_);
  }

  void step()
  {
    {
      const auto dt  = dt_;
      const auto dt2 = dt*Real{1.0/2.0};
      const auto dt3 = dt*Real{1.0/3.0};

#pragma omp parallel for schedule(runtime)
      for (int i = 0; i < nbody_; ++i)
      {
        posx_[i] += dt*(velx_[i] + dt2*(accx_[i] + dt3*jrkx_[i]));
        posy_[i] += dt*(vely_[i] + dt2*(accy_[i] + dt3*jrky_[i]));
        posz_[i] += dt*(velz_[i] + dt2*(accz_[i] + dt3*jrkz_[i]));

        velx_[i] += dt*(accx_[i] + dt2*jrkx_[i]);
        vely_[i] += dt*(accy_[i] + dt2*jrky_[i]);
        velz_[i] += dt*(accz_[i] + dt2*jrkz_[i]);
      }
    }

    ispc::compute_forces(
        nbody_,
        mass_.data(),
        posx_.data(),
        posy_.data(),
        posz_.data(),
        velx_.data(),
        vely_.data(),
        velz_.data(),
        accx1_.data(),
        accy1_.data(),
        accz1_.data(),
        jrkx1_.data(),
        jrky1_.data(),
        jrkz1_.data(),
        gpot_.data(),
        eps2_);

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
#pragma omp parallel for schedule(runtime) reduction(min:dt_min)
      for (int i = 0; i < nbody_; ++i)
      {
        /* interpolate snap & crackle */
        const auto Amx = accx1_[i] - accx_[i];
        const auto Amy = accy1_[i] - accy_[i];
        const auto Amz = accz1_[i] - accz_[i];

        const auto Jmx = h*(jrkx1_[i] - jrkx_[i]);
        const auto Jmy = h*(jrky1_[i] - jrky_[i]);
        const auto Jmz = h*(jrkz1_[i] - jrkz_[i]);

        const auto Jpx = h*(jrkx1_[i] + jrkx_[i]);
        const auto Jpy = h*(jrky1_[i] + jrky_[i]);
        const auto Jpz = h*(jrkz1_[i] + jrkz_[i]);

        auto snpx = f1 * Jmx;
        auto snpy = f1 * Jmy;
        auto snpz = f1 * Jmz;

        auto crkx = f2 * (Jpx - Amx);
        auto crky = f2 * (Jpy - Amy);
        auto crkz = f2 * (Jpz - Amz);

        snpx -= h*crkx;
        snpy -= h*crky;
        snpz -= h*crkz;

        /* correct */

        posx_[i] += dt4*snpx + dt5*crkx;
        posy_[i] += dt4*snpy + dt5*crky;
        posz_[i] += dt4*snpz + dt5*crkz;

        velx_[i] += dt3*snpx + dt4*crkx;
        vely_[i] += dt3*snpy + dt4*crky;
        velz_[i] += dt3*snpz + dt4*crkz;

        accx_[i] = accx1_[i];
        accy_[i] = accy1_[i];
        accz_[i] = accz1_[i];

        jrkx_[i] = jrkx1_[i];
        jrky_[i] = jrky1_[i];
        jrkz_[i] = jrkz1_[i];

        /* compute new timestep */

        const auto s0 = accx_[i]*accx_[i] + accy_[i]*accy_[i] + accz_[i]*accz_[i];
        const auto s1 = jrkx_[i]*jrkx_[i] + jrky_[i]*jrky_[i] + jrkz_[i]*jrkz_[i];
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
#pragma omp parallel for reduction(+:ekin,epot)
    for (int i = 0; i < nbody_; i++)
    {
      ekin += mass_[i] * (square(velx_[i]) + square(vely_[i]) + square(velz_[i])) * Real{0.5};
      epot += Real{0.5} * mass_[i] * gpot_[i];
    }

    return epot + ekin;
  }

};

void run(const int nbody, const Real t_end, const int dstep, const Real eta, const Real Q)
{
  using namespace std;
  Hermite4 h4(nbody,Q,eta);


  h4.step();

  const auto etot0 = h4.etot();
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

  auto param_pack = pack(argc, argv,
      param("number of bodies",     nbodies, "n", "nbodies"),
      param("integration time",     t_end,  "t", "tend"),
      param("print output every # steps" , dstep, "s", "dstep"),
      param("accuracy parameter",   eta,     "e", "eta"),
      param("virial ratio",         Q,       "Q", "")
      );

  cerr << param_pack.parse_all();

  assert(nbodies > 1);
  assert(t_end > 0);
  assert(eta > 0);
  assert(Q > 0);

  run(nbodies, t_end, dstep, eta, Q);

  return 0;
}

