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

#include "typeReal.h"


typedef Real<3> vec3;
struct Force
{
  vec3 acc, jrk;
  Real pot, null;
};

struct Predictor
{
  vec3 pos, vel;
};

static inline
void body_body_force(
    Force &fi,
    const Predictor &pi,
    const Predictor &pj,
    const Real mj,
    const Real eps2)
{
  const Real dx = pj.pos.x - pi.pos.x;
  const Real dy = pj.pos.y - pi.pos.y;
  const Real dz = pj.pos.z - pi.pos.z;

  const Real r2 =  dx*dx + dy*dy + dz*dz;
  const Real ds2 = r2 + eps2;

  const Real  inv_ds  = rsqrt(ds2);
  const Real  inv_ds2 = inv_ds*inv_ds;
  const Real minv_ds  = inv_ds  * mj;
  const Real minv_ds3 = inv_ds2 * minv_ds;

  fi.acc.x += minv_ds3 * dx;
  fi.acc.y += minv_ds3 * dy;
  fi.acc.z += minv_ds3 * dz;
  fi.pot   -= minv_ds * (Real)(r2 > 0.0);

  const Real dvx = pj.vel.x - pi.vel.x;
  const Real dvy = pj.vel.y - pi.vel.y;
  const Real dvz = pj.vel.z - pi.vel.z;
  const Real rv  = dx*dvx + dy*dvy + dz*dvz;

  const Real Jij = (Real)(-3.0) * (rv * inv_ds2 * minv_ds3);

  fi.jrk.x += minv_ds3*dvx + Jij*dx;
  fi.jrk.y += minv_ds3*dvy + Jij*dy;
  fi.jrk.z += minv_ds3*dvz + Jij*dz;
}

task void compute_forces_task(
    uniform const int     n,
    uniform const int nPerTask,
    uniform const Real mass[],
    uniform const Real posx[],
    uniform const Real posy[],
    uniform const Real posz[],
    uniform const Real velx[],
    uniform const Real vely[],
    uniform const Real velz[],
    uniform       Real accx[],
    uniform       Real accy[],
    uniform       Real accz[],
    uniform       Real jrkx[],
    uniform       Real jrky[],
    uniform       Real jrkz[],
    uniform       Real gpot[],
    const uniform Real eps2)
{
  const uniform int nibeg = taskIndex * nPerTask;
  const uniform int niend = min(n, nibeg + nPerTask);

  if (nibeg >= n)
    return;

  uniform Real shdata[7][programCount];

  assert((n%programCount) == 0);

  foreach (i = nibeg ... niend)
  {
    Force fi;
    fi.acc = (Real)0.0;
    fi.jrk = (Real)0.0;
    fi.pot = (Real)0.0;

    Predictor pi;
    pi.pos.x = posx[i];
    pi.pos.y = posy[i];
    pi.pos.z = posz[i];
    pi.vel.x = velx[i];
    pi.vel.y = vely[i];
    pi.vel.z = velz[i];

    for (uniform int jb = 0; jb < n; jb += programCount)
    {
      const int jp = jb + programIndex;
      shdata[0][programIndex] = posx[jp];
      shdata[1][programIndex] = posy[jp];
      shdata[2][programIndex] = posz[jp];
      shdata[3][programIndex] = mass[jp];
      shdata[4][programIndex] = velx[jp];
      shdata[5][programIndex] = vely[jp];
      shdata[6][programIndex] = velz[jp];

      for (uniform int j = 0; j < programCount; j++)
      {
        Predictor pj;
        pj.pos.x = shdata[0][j];
        pj.pos.y = shdata[1][j];
        pj.pos.z = shdata[2][j];
        pj.vel.x = shdata[4][j];
        pj.vel.y = shdata[5][j];
        pj.vel.z = shdata[6][j];
        const Real jmass  = shdata[3][j];
        body_body_force(fi,pi,pj,jmass,eps2);
      }
    }

    accx[i] = fi.acc.x;
    accy[i] = fi.acc.y;
    accz[i] = fi.acc.z;
    jrkx[i] = fi.jrk.x;
    jrky[i] = fi.jrk.y;
    jrkz[i] = fi.jrk.z;
    gpot[i] = fi.pot;
  }
}

export void compute_forces(
    uniform const int     n,
    uniform const Real mass[],
    uniform const Real posx[],
    uniform const Real posy[],
    uniform const Real posz[],
    uniform const Real velx[],
    uniform const Real vely[],
    uniform const Real velz[],
    uniform       Real accx[],
    uniform       Real accy[],
    uniform       Real accz[],
    uniform       Real jrkx[],
    uniform       Real jrky[],
    uniform       Real jrkz[],
    uniform       Real gpot[],
    const uniform Real eps2)
{
  const uniform int nPerTask = min(128,programCount*8);
  const uniform int nTask = (n+nPerTask-1)/nPerTask;

  launch [nTask]  compute_forces_task(
      n, nPerTask,
      mass,
      posx,posy,posz,
      velx,vely,velz,
      accx,accy,accz,
      jrkx,jrky,jrkz,
      gpot,eps2);
}
