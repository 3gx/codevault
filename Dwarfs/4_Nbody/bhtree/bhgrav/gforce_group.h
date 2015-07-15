#ifndef __GFORCE_GROUP_H__
#define __GFORCE_GROUP_H__

#if 0
#define DRY_RUN
#endif

#if 1
#define QUADRUPOLE
#endif

inline bool split_node(
    const float4 &cellCoM,     /* contains (x,y,z,w), where w is an opening criteria distance */
    const float4 &groupCentre,
    const float4 &groupSize) const

{
  const float4  dr = fabs(groupCentre - cellCoM) - groupSize;
  const float4 adr = dr + fabs(dr);  /* 0 if dr < 0, otherwise 2*dr */
  const float  ds2 = adr.norm2()*0.25f;

  return ds2 <= cellCoM.w();
}

  template<const int Ng>
std::pair<unsigned int, unsigned int> gForce(const GroupT<Ng> &group, float4 force[Ng]) const
{
  const int Ncell = NLEAF*16;
  const int Nptcl = NLEAF*64;
  float4      ptcl_list[Nptcl*2] __attribute__ ((aligned(32)));
  fMultipole  cell_list[Ncell*2] __attribute__ ((aligned(32)));
  int nc = 0, np = 0;
  int np_tot = 0, nc_tot = 0;

  gForce<true, Ncell, Nptcl>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot);

  return std::make_pair(np_tot*group.nb(), nc_tot*group.nb());
}

template<const bool ROOT, const int Nc, const int Np, const int Ng>
void gForce(
    const GroupT<Ng> &group,
    float4         force[Ng  ],         /* maximal number of particles in a group  */
    fMultipole cell_list[Nc*2], int &nc, /* list for particle-cell interactions     */
    float4     ptcl_list[Np*2], int &np, /* list for particle-particle interactions */
    int &np_tot, int &nc_tot,
    const  float4 &groupCentre = 0.0, const float4 &groupSize = 0.0, const int addr = 0) const
{
  if (ROOT)
  {
    assert(isTreeReady());
    assert(Np >= 8*NLEAF);
    assert(Nc >  8);
    const boundary bnd = group.outerBoundary();
    const float4 groupCentre =  bnd.center().           to_float4();
    const float4 groupSize   = (bnd.hlen  ()*(real)2.0).to_float4();
    nc = np = 0;
    nc_tot = np_tot = 0;
    for (int i = 0; i < Ng; i++)
      force[i] = 0.0f;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
        if (split_node(cellCoM[cellList[k].id()], groupCentre, groupSize))
          gForce<false, Nc, Np>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot,
              groupCentre, groupSize, k);

    np = particle_particle       <Np>(group, force, ptcl_list, np);
    np = particle_particle_scalar<Np>(group, force, ptcl_list, np);
    nc = particle_cell           <Nc>(group, force, cell_list, nc);
    nc = particle_cell_scalar    <Nc>(group, force, cell_list, nc);
  }
  else
  {
    const Cell cell = cellList[addr];
    if (split_node(cellCoM[cell.id()], groupCentre, groupSize))
    {
      if (cell.isNode())
      {
        for (int k = 0; k < 8; k++)
          if (!cellList[cell.addr()+k].isEmpty())
            gForce<false, Nc, Np>(group, force, cell_list, nc, ptcl_list, np, np_tot, nc_tot,
                groupCentre, groupSize, cell.addr()+k);
      }
      else
      {
        const Leaf &leaf = leafList[cell.leafIdx()];
        for (int i = 0; i < leaf.nb(); i++)
          ptcl_list[np++] = leaf[i].pos_mass();
        assert(np <= 2*Np);
        np_tot += leaf.nb();
      }
    }
    else
    {
      cell_list[nc++] = multipoleList[cell.id()];
      nc_tot++;
      assert(nc <= 2*Nc);
    }

    if (np >= Np) np = particle_particle<Np>(group, force, ptcl_list, np);
    if (nc >= Nc) nc = particle_cell    <Nc>(group, force, cell_list, nc);
  }
}

template<const int Np, const int Ng>
int particle_particle(
    const GroupT<Ng> &igroup, 
    float4     force[Ng  ],
    float4 ptcl_list[Np*2], const int np) const
{
#ifdef DRY_RUN
  return 0;
#endif
#ifdef __AVX_H__
#if 0
  {
    const GroupT<Ng> group = igroup;
    const int ni = group.nb();
    const int nj = np & (-8);

    const _v8sf *ib = (const _v8sf*)&group[0];
    const _v4sf *jb = (const _v4sf*)&ptcl_list[np - nj];

    const _v8sf veps2 = v8sf(eps2);
    for (int i = 0; i < ni; i++)
    {
      const _v8sf ip  = *(ib + i);
      const _v8sf ipx = __bcast<0>(ip);
      const _v8sf ipy = __bcast<1>(ip);
      const _v8sf ipz = __bcast<2>(ip);

      _v8sf fx   = v8sf(0.0f);
      _v8sf fy   = v8sf(0.0f);
      _v8sf fz   = v8sf(0.0f);
      _v8sf gpot = v8sf(0.0f);
      assert((nj&7) == 0);
      for (int j = 0; j < nj; j += 8)
      {
        const _v8sf jp04 = pack2ymm(*(jb+j+0), *(jb+j+4));
        const _v8sf jp15 = pack2ymm(*(jb+j+1), *(jb+j+5));
        const _v8sf jp26 = pack2ymm(*(jb+j+2), *(jb+j+6));
        const _v8sf jp37 = pack2ymm(*(jb+j+3), *(jb+j+7));

        const _v8sf xy02 = __builtin_ia32_unpcklps256(jp04, jp26);
        const _v8sf xy13 = __builtin_ia32_unpcklps256(jp15, jp37);
        const _v8sf zw02 = __builtin_ia32_unpckhps256(jp04, jp26);
        const _v8sf zw13 = __builtin_ia32_unpckhps256(jp15, jp37);
        const _v8sf  jpx = __builtin_ia32_unpcklps256(xy02, xy13);
        const _v8sf  jpy = __builtin_ia32_unpckhps256(xy02, xy13);
        const _v8sf  jpz = __builtin_ia32_unpcklps256(zw02, zw13);
        const _v8sf  jpm = __builtin_ia32_unpckhps256(zw02, zw13);

        const _v8sf dx = jpx - ipx;
        const _v8sf dy = jpy - ipy;
        const _v8sf dz = jpz - ipz;
        const _v8sf r2 = dx*dx + dy*dy + dz*dz + veps2;

        const _v8sf  rinv  = __builtin_ia32_rsqrtps256(r2);
        const _v8sf mrinv  = rinv  *  jpm;
        const _v8sf  rinv2 = rinv  *  rinv;
        const _v8sf mrinv3 = rinv2 * mrinv;

        fx   += mrinv3 * dx;
        fy   += mrinv3 * dy;
        fz   += mrinv3 * dz;
        gpot += mrinv;
      }
      gpot = -gpot;

      const _v4sf f4x = __reduce_v8sf(fx);
      const _v4sf f4y = __reduce_v8sf(fz);
      const _v4sf f4z = __reduce_v8sf(fy);
      const _v4sf f4w = __reduce_v8sf(gpot);

      force[i] = force[i] + float4(
          __builtin_ia32_vec_ext_v4sf(f4x, 0),
          __builtin_ia32_vec_ext_v4sf(f4y, 1),
          __builtin_ia32_vec_ext_v4sf(f4z, 2),
          __builtin_ia32_vec_ext_v4sf(f4w, 3));
    }
    return np - nj;
  }
#else
  {
    const GroupT<Ng> group = igroup;
    const _v8sf *ib = (const _v8sf*)&group[0];
    const _v4sf *jb = (const _v4sf*)ptcl_list;

    const int ni = group.nb();
    const int nj = np;

    const _v8sf veps2 = v8sf(eps2);
    for (int i = 0; i < ni; i += 8)
    {
      const _v8sf ip04 = __mergelo(*(ib+i+0), *(ib+i+4));
      const _v8sf ip15 = __mergelo(*(ib+i+1), *(ib+i+5));
      const _v8sf ip26 = __mergelo(*(ib+i+2), *(ib+i+6));
      const _v8sf ip37 = __mergelo(*(ib+i+3), *(ib+i+7));

      const _v8sf xy02 = __builtin_ia32_unpcklps256(ip04, ip26);
      const _v8sf xy13 = __builtin_ia32_unpcklps256(ip15, ip37);
      const _v8sf zw02 = __builtin_ia32_unpckhps256(ip04, ip26);
      const _v8sf zw13 = __builtin_ia32_unpckhps256(ip15, ip37);
      const _v8sf  ipx = __builtin_ia32_unpcklps256(xy02, xy13);
      const _v8sf  ipy = __builtin_ia32_unpckhps256(xy02, xy13);
      const _v8sf  ipz = __builtin_ia32_unpcklps256(zw02, zw13);

      _v8sf fx   = v8sf(0.0f);
      _v8sf fy   = v8sf(0.0f);
      _v8sf fz   = v8sf(0.0f);
      _v8sf gpot = v8sf(0.0f);
      for (int j = 0; j < nj; j++)
      {
        const _v8sf jp  = __bcast0(jb+j);
        const _v8sf jpx = __builtin_ia32_shufps256(jp, jp, 0x00);
        const _v8sf jpy = __builtin_ia32_shufps256(jp, jp, 0x55);
        const _v8sf jpz = __builtin_ia32_shufps256(jp, jp, 0xAA);
        const _v8sf jpm = __builtin_ia32_shufps256(jp, jp, 0xFF);

        const _v8sf dx = jpx - ipx;
        const _v8sf dy = jpy - ipy;
        const _v8sf dz = jpz - ipz;
        const _v8sf r2 = dx*dx + dy*dy + dz*dz + veps2;

        const _v8sf  rinv  = __builtin_ia32_rsqrtps256(r2);
        const _v8sf mrinv  = rinv  *  jpm;
        const _v8sf  rinv2 = rinv  *  rinv;
        const _v8sf mrinv3 = rinv2 * mrinv;

        fx   += mrinv3 * dx;
        fy   += mrinv3 * dy;
        fz   += mrinv3 * dz;
        gpot += mrinv;
      }
      gpot = -gpot;

      const _v8sf t0  = __builtin_ia32_unpcklps256(fx,   fz);
      const _v8sf t1  = __builtin_ia32_unpcklps256(fy, gpot);
      const _v8sf t2  = __builtin_ia32_unpckhps256(fx,   fz);
      const _v8sf t3  = __builtin_ia32_unpckhps256(fy, gpot);
      const _v8sf f04 = __builtin_ia32_unpcklps256(t0, t1);
      const _v8sf f15 = __builtin_ia32_unpckhps256(t0, t1);
      const _v8sf f26 = __builtin_ia32_unpcklps256(t2, t3);
      const _v8sf f37 = __builtin_ia32_unpckhps256(t2, t3);

#if 0
      _v8sf* vforce = (_v8sf*)&force[i];
      *(vforce + 0) += __merge<0,0>(f04, f15);
      *(vforce + 1) += __merge<0,0>(f26, f37);
      *(vforce + 2) += __merge<1,1>(f04, f15);
      *(vforce + 3) += __merge<1,1>(f26, f37);
#else
      _v4sf* vforce = (_v4sf*)&force[i];
      *(vforce + 0) += __extract<0>(f04);
      *(vforce + 4) += __extract<1>(f04);
      *(vforce + 1) += __extract<0>(f15);
      *(vforce + 5) += __extract<1>(f15);
      *(vforce + 2) += __extract<0>(f26);
      *(vforce + 6) += __extract<1>(f26);
      *(vforce + 3) += __extract<0>(f37);
      *(vforce + 7) += __extract<1>(f37);
#endif
    }
    return 0;
  }
#endif
#elif defined __SSE_H__
  {
    const GroupT<Ng> &group = igroup;
    const _v4sf *ib = (const _v4sf*)&group[0];
    const _v4sf *jb = (const _v4sf*)ptcl_list;

    const int ni = group.nb();
    const int nj = np;

    const _v4sf veps2 = v4sf(eps2);
    for (int i = 0; i < ni; i += 4)
    {
      const int i2 = i<<1;
      const _v4sf ip0 = *(ib + i2 + 0);
      const _v4sf ip1 = *(ib + i2 + 2);
      const _v4sf ip2 = *(ib + i2 + 4);
      const _v4sf ip3 = *(ib + i2 + 6);

      const _v4sf t0  = __builtin_ia32_unpcklps(ip0, ip2);
      const _v4sf t1  = __builtin_ia32_unpckhps(ip0, ip2);
      const _v4sf t2  = __builtin_ia32_unpcklps(ip1, ip3);
      const _v4sf t3  = __builtin_ia32_unpckhps(ip1, ip3);
      const _v4sf ipx = __builtin_ia32_unpcklps(t0,  t2);
      const _v4sf ipy = __builtin_ia32_unpckhps(t0,  t2);
      const _v4sf ipz = __builtin_ia32_unpcklps(t1,  t3);

      _v4sf fx   = v4sf(0.0f);
      _v4sf fy   = v4sf(0.0f);
      _v4sf fz   = v4sf(0.0f);
      _v4sf gpot = v4sf(0.0f);
      for (int j = 0; j < nj; j++)
      {
        const _v4sf jp = *(jb + j);
        const _v4sf jpx = __builtin_ia32_shufps(jp, jp, 0x00);
        const _v4sf jpy = __builtin_ia32_shufps(jp, jp, 0x55);
        const _v4sf jpz = __builtin_ia32_shufps(jp, jp, 0xAA);
        const _v4sf jpm = __builtin_ia32_shufps(jp, jp, 0xFF);

        const _v4sf dx = jpx - ipx;
        const _v4sf dy = jpy - ipy;
        const _v4sf dz = jpz - ipz;
        const _v4sf r2 = dx*dx + dy*dy + dz*dz + veps2;

        const _v4sf  rinv  = __builtin_ia32_rsqrtps(r2);
        const _v4sf mrinv  = rinv  *  jpm;
        const _v4sf  rinv2 = rinv  *  rinv;
        const _v4sf mrinv3 = rinv2 * mrinv;

        fx   += mrinv3 * dx;
        fy   += mrinv3 * dy;
        fz   += mrinv3 * dz;
        gpot += mrinv;
      }
      gpot = -gpot;

      {
        const _v4sf t0 = __builtin_ia32_unpcklps(fx,   fz);
        const _v4sf t1 = __builtin_ia32_unpcklps(fy, gpot);
        const _v4sf t2 = __builtin_ia32_unpckhps(fx,   fz);
        const _v4sf t3 = __builtin_ia32_unpckhps(fy, gpot);
        const _v4sf f0 = __builtin_ia32_unpcklps(t0, t1);
        const _v4sf f1 = __builtin_ia32_unpckhps(t0, t1);
        const _v4sf f2 = __builtin_ia32_unpcklps(t2, t3);
        const _v4sf f3 = __builtin_ia32_unpckhps(t2, t3);

        _v4sf *vforce = (_v4sf*)&force[i];
        *(vforce + 0) += f0;
        *(vforce + 1) += f1;
        *(vforce + 2) += f2;
        *(vforce + 3) += f3;
      }
    }
    return 0;
  }
#else  /* !__AVX_H__ */
  return particle_particle_scalar<Np>(igroup, force, ptcl_list, np);
#endif
}

template<const int Nc, const int Ng>
int particle_cell(
    const GroupT<Ng> &igroup,
    float4         force[Ng  ],       
    fMultipole cell_list[Nc*2], const int nc) const
{
#ifdef DRY_RUN
  return 0;
#endif
#ifdef __AVX_H__
  const GroupT<Ng> group = igroup;
  const _v8sf *ib = (const _v8sf*)&group[0];
  const _v4sf *jb = (const _v4sf*)cell_list;

  const int ni  = group.nb();
  const int nj  = nc;
  const int nj3 = nj * 3;

  const _v8sf veps2 = v8sf(eps2);
  for (int i = 0; i < ni; i += 8)
  {
    const _v8sf ip04 = __mergelo(*(ib+i+0), *(ib+i+4));
    const _v8sf ip15 = __mergelo(*(ib+i+1), *(ib+i+5));
    const _v8sf ip26 = __mergelo(*(ib+i+2), *(ib+i+6));
    const _v8sf ip37 = __mergelo(*(ib+i+3), *(ib+i+7));

    const _v8sf xy02 = __builtin_ia32_unpcklps256(ip04, ip26);
    const _v8sf xy13 = __builtin_ia32_unpcklps256(ip15, ip37);
    const _v8sf zw02 = __builtin_ia32_unpckhps256(ip04, ip26);
    const _v8sf zw13 = __builtin_ia32_unpckhps256(ip15, ip37);
    const _v8sf  ipx = __builtin_ia32_unpcklps256(xy02, xy13);
    const _v8sf  ipy = __builtin_ia32_unpckhps256(xy02, xy13);
    const _v8sf  ipz = __builtin_ia32_unpcklps256(zw02, zw13);

    _v8sf fx   = v8sf(0.0f);
    _v8sf fy   = v8sf(0.0f);
    _v8sf fz   = v8sf(0.0f);
    _v8sf gpot = v8sf(0.0f);
    for (int j = 0; j < nj3; j += 3)
    {
      const _v8sf jp  = __bcast0(jb+j);
      const _v8sf jpx = __builtin_ia32_shufps256(jp, jp, 0x00);
      const _v8sf jpy = __builtin_ia32_shufps256(jp, jp, 0x55);
      const _v8sf jpz = __builtin_ia32_shufps256(jp, jp, 0xAA);
      const _v8sf jpm = __builtin_ia32_shufps256(jp, jp, 0xFF);

      const _v8sf dx = jpx - ipx;
      const _v8sf dy = jpy - ipy;
      const _v8sf dz = jpz - ipz;
      const _v8sf r2 = dx*dx + dy*dy + dz*dz + veps2;

      const _v8sf  rinv  = __builtin_ia32_rsqrtps256(r2);
      const _v8sf mrinv  = rinv  *  jpm;
      const _v8sf  rinv2 = rinv  *  rinv;
      const _v8sf mrinv3 = rinv2 * mrinv;

      fx   += mrinv3 * dx;
      fy   += mrinv3 * dy;
      fz   += mrinv3 * dz;
      gpot += mrinv;

      /* 21 flops */

#ifdef QUADRUPOLE
      const _v8sf rinv4 = rinv2*rinv2;
      const _v8sf rinv5 = rinv4*rinv;
      const _v8sf rinv7 = rinv5*rinv2;

      const _v8sf Q1  = __bcast0(jb+j+1);
      const _v8sf Qxx = __builtin_ia32_shufps256(Q1, Q1, 0x00);
      const _v8sf Qyy = __builtin_ia32_shufps256(Q1, Q1, 0x55);
      const _v8sf Qzz = __builtin_ia32_shufps256(Q1, Q1, 0xAA);

      const _v8sf Q2  = __bcast0(jb+j+2);
      const _v8sf Qxy = __builtin_ia32_shufps256(Q2, Q2, 0x00);
      const _v8sf Qxz = __builtin_ia32_shufps256(Q2, Q2, 0x55);
      const _v8sf Qyz = __builtin_ia32_shufps256(Q2, Q2, 0xAA);

      const _v8sf Qrx = Qxx*dx + Qxy*dy + Qxz*dz;
      const _v8sf Qry = Qxy*dx + Qyy*dy + Qyz*dz;
      const _v8sf Qrz = Qxz*dx + Qyz*dy + Qzz*dz;

      const _v8sf rQr = Qrx*dx + Qry*dy + Qrz*dz;

      const _v8sf c1 = v8sf(2.5f)*rinv7*rQr;
      fx   += c1*dx;
      fy   += c1*dy;
      fz   += c1*dz;
      gpot += v8sf(0.5f)*rinv5*rQr;

      const _v8sf c2 = -rinv5;
      fx += c2*Qrx;
      fy += c2*Qry;
      fz += c2*Qrz;

      /* 40 flops */
#endif
    }
    gpot = -gpot;

    const _v8sf t0  = __builtin_ia32_unpcklps256(fx,   fz);
    const _v8sf t1  = __builtin_ia32_unpcklps256(fy, gpot);
    const _v8sf t2  = __builtin_ia32_unpckhps256(fx,   fz);
    const _v8sf t3  = __builtin_ia32_unpckhps256(fy, gpot);
    const _v8sf f04 = __builtin_ia32_unpcklps256(t0, t1);
    const _v8sf f15 = __builtin_ia32_unpckhps256(t0, t1);
    const _v8sf f26 = __builtin_ia32_unpcklps256(t2, t3);
    const _v8sf f37 = __builtin_ia32_unpckhps256(t2, t3);

#if 0
    _v8sf* vforce = (_v8sf*)&force[i];
    *(vforce + 0) += __merge<0,0>(f04, f15);
    *(vforce + 1) += __merge<0,0>(f26, f37);
    *(vforce + 2) += __merge<1,1>(f04, f15);
    *(vforce + 3) += __merge<1,1>(f26, f37);
#else
    _v4sf* vforce = (_v4sf*)&force[i];
    *(vforce + 0) += __extract<0>(f04);
    *(vforce + 4) += __extract<1>(f04);
    *(vforce + 1) += __extract<0>(f15);
    *(vforce + 5) += __extract<1>(f15);
    *(vforce + 2) += __extract<0>(f26);
    *(vforce + 6) += __extract<1>(f26);
    *(vforce + 3) += __extract<0>(f37);
    *(vforce + 7) += __extract<1>(f37);
#endif
  }
  return 0;
#elif defined __SSE_H__
  {
    const GroupT<Ng> &group = igroup;
    const _v4sf *ib = (const _v4sf*)&group[0];
    const _v4sf *jb = (const _v4sf*)cell_list;

    const int ni  = group.nb();
    const int nj  = nc;
    const int nj3 = nj * 3;

    const _v4sf veps2 = v4sf(eps2);
    for (int i = 0; i < ni; i += 4)
    {
      const int i2 = i<<1;
      const _v4sf ip0 = *(ib + i2 + 0);
      const _v4sf ip1 = *(ib + i2 + 2);
      const _v4sf ip2 = *(ib + i2 + 4);
      const _v4sf ip3 = *(ib + i2 + 6);

      const _v4sf t0  = __builtin_ia32_unpcklps(ip0, ip2);
      const _v4sf t1  = __builtin_ia32_unpckhps(ip0, ip2);
      const _v4sf t2  = __builtin_ia32_unpcklps(ip1, ip3);
      const _v4sf t3  = __builtin_ia32_unpckhps(ip1, ip3);
      const _v4sf ipx = __builtin_ia32_unpcklps(t0,  t2);
      const _v4sf ipy = __builtin_ia32_unpckhps(t0,  t2);
      const _v4sf ipz = __builtin_ia32_unpcklps(t1,  t3);

      _v4sf fx   = v4sf(0.0f);
      _v4sf fy   = v4sf(0.0f);
      _v4sf fz   = v4sf(0.0f);
      _v4sf gpot = v4sf(0.0f);
      for (int j = 0; j < nj3; j += 3)
      {
        const _v4sf jp = *(jb + j);
        const _v4sf jpx = __builtin_ia32_shufps(jp, jp, 0x00);
        const _v4sf jpy = __builtin_ia32_shufps(jp, jp, 0x55);
        const _v4sf jpz = __builtin_ia32_shufps(jp, jp, 0xAA);
        const _v4sf jpm = __builtin_ia32_shufps(jp, jp, 0xFF);

        const _v4sf dx = jpx - ipx;
        const _v4sf dy = jpy - ipy;
        const _v4sf dz = jpz - ipz;
        const _v4sf r2 = dx*dx + dy*dy + dz*dz + veps2;

        const _v4sf  rinv  = __builtin_ia32_rsqrtps(r2);
        const _v4sf mrinv  = rinv  *  jpm;
        const _v4sf  rinv2 = rinv  *  rinv;
        const _v4sf mrinv3 = rinv2 * mrinv;

        fx   += mrinv3 * dx;
        fy   += mrinv3 * dy;
        fz   += mrinv3 * dz;
        gpot += mrinv;

#ifdef QUADRUPOLE
        const _v4sf rinv4 = rinv2*rinv2;
        const _v4sf rinv5 = rinv4*rinv;
        const _v4sf rinv7 = rinv5*rinv2;

        const _v4sf Q1  = *(jb+j+1);
        const _v4sf Qxx = __builtin_ia32_shufps(Q1, Q1, 0x00);
        const _v4sf Qyy = __builtin_ia32_shufps(Q1, Q1, 0x55);
        const _v4sf Qzz = __builtin_ia32_shufps(Q1, Q1, 0xAA);
        const _v4sf Q2  = *(jb+j+2);
        const _v4sf Qxy = __builtin_ia32_shufps(Q2, Q2, 0x00);
        const _v4sf Qxz = __builtin_ia32_shufps(Q2, Q2, 0x55);
        const _v4sf Qyz = __builtin_ia32_shufps(Q2, Q2, 0xAA);

        const _v4sf Qrx = Qxx*dx + Qxy*dy + Qxz*dz;
        const _v4sf Qry = Qxy*dx + Qyy*dy + Qyz*dz;
        const _v4sf Qrz = Qxz*dx + Qyz*dy + Qzz*dz;

        const _v4sf rQr = Qrx*dx + Qry*dy + Qrz*dz;

        const _v4sf c1 = v4sf(2.5f)*rinv7*rQr;
        fx   += c1*dx;
        fy   += c1*dy;
        fz   += c1*dz;
        gpot += v4sf(0.5f)*rinv5*rQr;

        const _v4sf c2 = -rinv5;
        fx += c2*Qrx;
        fy += c2*Qry;
        fz += c2*Qrz;

        /* 40 flops */
#endif
      }
      gpot = -gpot;

      {
        const _v4sf t0 = __builtin_ia32_unpcklps(fx,   fz);
        const _v4sf t1 = __builtin_ia32_unpcklps(fy, gpot);
        const _v4sf t2 = __builtin_ia32_unpckhps(fx,   fz);
        const _v4sf t3 = __builtin_ia32_unpckhps(fy, gpot);
        const _v4sf f0 = __builtin_ia32_unpcklps(t0, t1);
        const _v4sf f1 = __builtin_ia32_unpckhps(t0, t1);
        const _v4sf f2 = __builtin_ia32_unpcklps(t2, t3);
        const _v4sf f3 = __builtin_ia32_unpckhps(t2, t3);

        _v4sf *vforce = (_v4sf*)&force[i];
        *(vforce + 0) += f0;
        *(vforce + 1) += f1;
        *(vforce + 2) += f2;
        *(vforce + 3) += f3;
      }
    }
    return 0;
  }
#else  /* __AVX_H__ */
  return particle_cell_scalar<Nc>(igroup, force, cell_list, nc);
#endif
}

/**********************/
/*** SCALAR VERSION ***/
/**********************/

template<const int Np, const int Ng>
int particle_particle_scalar(
    const GroupT<Ng> &group, 
    float4     force[Ng  ],
    float4 ptcl_list[Np*2], const int np) const
{
  if (np == 0) return 0;
  const int ni = group.nb();
  const int nj = np;
  for (int i = 0; i < ni; i++)
  {
    const float4 ip = group[i].pos_h();
    for (int j = 0; j < nj; j++)
    {
      const float4 jp = ptcl_list[j];
      const float4 dr = jp - ip;
      const float  r2 = dr.norm2() + eps2;
      const float  mj = jp.w();

      const float  rinv  = 1.0f/std::sqrt(r2);
      const float mrinv  = rinv*mj;
      const float mrinv3 = rinv*rinv*mrinv;

      float4 acc = dr * float4(mrinv3);
      acc.w() = -mrinv;

      force[i] = force[i] + acc;
    }
  }
  return 0;
}

template<const int Nc, const int Ng>
int particle_cell_scalar(
    const GroupT<Ng> &group,
    float4         force[Ng  ],       
    fMultipole cell_list[Nc*2], const int nc) const
{
  if (nc == 0) return 0;
  const int ni = group.nb();
  const int nj = nc;

  for (int i = 0; i < ni; i++)
  {
    const float4 ip = group[i].pos_mass();
    for (int j = 0; j < nj; j++)
    {
      const fMultipole &multipole = cell_list[j];
      const fMonopole   &m =   multipole.monopole();

      const float4 jp = *(float4*)&m;
      const float4 dr = jp - ip;
      const float  r2 = dr.norm2() + eps2;
      assert(r2 > 0.0f);
      const float  mj = jp.w();

      const float  rinv  = 1.0f/std::sqrt(r2);
      const float mrinv  = rinv*mj;
      const float mrinv3 = rinv*rinv*mrinv;

      float4 acc = dr * float4(mrinv3);
      acc.w() = -mrinv;

#ifdef QUADRUPOLE
      const fQuadrupole &q = multipole.quadrupole();
      const float rinv2 = rinv*rinv;
      const float rinv4 = rinv2*rinv2;
      const float rinv5 = rinv4*rinv;
      const float rinv7 = rinv5*rinv2;

      const float Qrx = (float4(q.xx(), q.xy(), q.xz(), 0.0f)*dr).reduce();
      const float Qry = (float4(q.xy(), q.yy(), q.yz(), 0.0f)*dr).reduce();
      const float Qrz = (float4(q.xz(), q.yz(), q.zz(), 0.0f)*dr).reduce();

      const float rQr = (float4(Qrx, Qry, Qrz, 0.0f)*dr).reduce();

      float4 qacc1 = float4(2.5f*rinv7*rQr)*dr;
      qacc1.w() = (-0.5f)*rinv5*rQr;

      qacc1 += float4(-rinv5)*float4(Qrx, Qry, Qrz, 0.0f);
      acc = acc +  qacc1;
#endif

      force[i] += acc;
    }
  }
  return 0;
}



#endif /* __GFORCE_GROUP_H__ */
