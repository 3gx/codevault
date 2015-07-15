#ifndef __RANGE_SEARCH_GROUP_H__
#define __RANGE_SEARCH_GROUP_H__

template<const bool ROOT, const int N>  /* must be ROOT = true on the root node (first call) */
void range_search(
    int nb[N],
    const GroupT<N> &igroup, const boundary &ibnd = boundary(), const int addr = 0) const
{
  if (ROOT)
  {
    assert(isTreeReady());
    const boundary ibnd = igroup.outerBoundary();
    for (int i = 0; i < N; i++)
      nb[i] = 0;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
        if (!not_overlapped(ibnd, bndsList[cellList[k].id()].inner()))
          range_search<false>(nb, igroup, ibnd, k);
  }
  else
  {
    const Cell cell = cellList[addr];
    if (cell.isNode())
    {
      for (int k = 0; k < 8; k++)
        if (!cellList[cell.addr()+k].isEmpty())
          if (!not_overlapped(ibnd, bndsList[cellList[cell.addr()+k].id()].inner()))
            range_search<false>(nb, igroup, ibnd, cell.addr()+k);
    }
    else
    {
      const Leaf      leaf  = leafList[cell.leafIdx()];
      const boundary &leafBnd = bndsList[cell.     id()].inner();
#ifdef __AVX_H__
      const GroupT<N> group = igroup;
      const _v4sf  jmin = leafBnd.min;
      const _v4sf  jmax = leafBnd.max;
      const int   ni = igroup.nb();
      const int   nj =   leaf.nb();
      const _v8sf *ib = (const _v8sf*)& group[0];
      const _v8sf *jb = (const _v8sf*)&  leaf[0];
      for (int i = 0; i < ni; i += 8)
      {
        asm("#eg01");
        const _v8sf ip04 = __mergelo(*(ib+i+0), *(ib+i+4));
        const _v8sf ip15 = __mergelo(*(ib+i+1), *(ib+i+5));
        const _v8sf ip26 = __mergelo(*(ib+i+2), *(ib+i+6));
        const _v8sf ip37 = __mergelo(*(ib+i+3), *(ib+i+7));

        /* check if these i-particles overlap with the leaf */

        const _v8sf  h04 = __bcast2<3>(ip04);
        const _v8sf  h15 = __bcast2<3>(ip15);
        const _v8sf  h26 = __bcast2<3>(ip26);
        const _v8sf  h37 = __bcast2<3>(ip37);
        const _v8sf imint = __builtin_ia32_minps256(
            __builtin_ia32_minps256(ip04-h04, ip15-h15),
            __builtin_ia32_minps256(ip26-h26, ip37-h37));
        const _v8sf imaxt = __builtin_ia32_maxps256(
            __builtin_ia32_maxps256(ip04+h04, ip15+h15),
            __builtin_ia32_maxps256(ip26+h26, ip37+h37));

        const _v4sf imin = __builtin_ia32_minps(
            __extract<0>(imint), __extract<1>(imint));
        const _v4sf imax = __builtin_ia32_maxps(
            __extract<0>(imaxt), __extract<1>(imaxt));

        const bool skip    = 
          __builtin_ia32_movmskps(__builtin_ia32_orps(
                __builtin_ia32_cmpltps(jmax, imin),
                __builtin_ia32_cmpltps(imax, jmin))) & 7;
        if (skip && i+7 < ni) continue;

        /* they do overlap, now proceed to the interaction part */

        const _v8sf xy02 = __builtin_ia32_unpcklps256(ip04, ip26);
        const _v8sf xy13 = __builtin_ia32_unpcklps256(ip15, ip37);
        const _v8sf zw02 = __builtin_ia32_unpckhps256(ip04, ip26);
        const _v8sf zw13 = __builtin_ia32_unpckhps256(ip15, ip37);
        const _v8sf xxxx = __builtin_ia32_unpcklps256(xy02, xy13);
        const _v8sf yyyy = __builtin_ia32_unpckhps256(xy02, xy13);
        const _v8sf zzzz = __builtin_ia32_unpcklps256(zw02, zw13);
        const _v8sf wwww = __builtin_ia32_unpckhps256(zw02, zw13);
        const _v8sf ipx = xxxx;
        const _v8sf ipy = yyyy;
        const _v8sf ipz = zzzz;
        const _v8sf iph = wwww;
        const _v8sf iph2 = iph*iph;

        _v8sf fnb = v8sf(0.0f);
        for (int j = 0; j < nj; j++)
        {
          const _v8sf jp = *(jb + j);

          const _v8sf jpx = __bcast<0>(jp);
          const _v8sf jpy = __bcast<1>(jp);
          const _v8sf jpz = __bcast<2>(jp);

          const _v8sf dx = jpx - ipx;
          const _v8sf dy = jpy - ipy;
          const _v8sf dz = jpz - ipz;
          const _v8sf r2 = dx*dx + dy*dy + dz*dz;

          const _v8sf mask = __builtin_ia32_cmpps256(r2, iph2, 1);
#if 0
          const int imask = __builtin_ia32_movmskps256(mask);
          if (imask == 0) continue;
#endif
#if 0 
          fnb += v8sf(1.0f);   /* this is to count number of ptcl processed */
#else
          fnb += __builtin_ia32_andps256(v8sf(1.0f), mask);
#endif
        }
        const _v8si inb = __builtin_ia32_cvtps2dq256(fnb);
        *(_v8si*)&nb[i] += inb;
        asm("#eg02");
      }
#elif defined __SSE_H__
      const _v4sf  jmin = leafBnd.min;
      const _v4sf  jmax = leafBnd.max;
      const int   ni = igroup.nb();
      const int   nj =   leaf.nb();
      const _v4sf *ib = (const _v4sf*)&igroup[0];
      const _v4sf *jb = (const _v4sf*)&  leaf[0];
      const int nj2 = nj<<1;
      for (int i = 0; i < ni; i += 4)
      {
        asm("#eg01");
        const int i2 = i<<1;
        const _v4sf ip0 = *(ib + i2 + 0);
        const _v4sf ip1 = *(ib + i2 + 2);
        const _v4sf ip2 = *(ib + i2 + 4);
        const _v4sf ip3 = *(ib + i2 + 6);

        /* check if these i-particles overlap with the leaf */

        const _v4sf h0  = __builtin_ia32_shufps(ip0, ip0, 0xFF);
        const _v4sf h1  = __builtin_ia32_shufps(ip1, ip1, 0xFF);
        const _v4sf h2  = __builtin_ia32_shufps(ip2, ip2, 0xFF);
        const _v4sf h3  = __builtin_ia32_shufps(ip3, ip3, 0xFF);

        /*   0     1     2     3   */
        /*  0x00 ,0x55, 0xaa, 0xff */

        const _v4sf imin    = __builtin_ia32_minps(__builtin_ia32_minps(ip0-h0, ip1-h1), __builtin_ia32_minps(ip2-h2,ip3-h3));
        const _v4sf imax    = __builtin_ia32_maxps(__builtin_ia32_maxps(ip0+h0, ip1+h1), __builtin_ia32_maxps(ip2+h2,ip3+h3));

        const bool skip    = __builtin_ia32_movmskps(__builtin_ia32_orps(
              __builtin_ia32_cmpltps(jmax, imin),
              __builtin_ia32_cmpltps(imax, jmin))) & 7;
        if (skip && i+4 < ni) continue;

        /* they do overlap, now proceed to the interaction part */

        const _v4sf t0 = __builtin_ia32_unpcklps(ip0, ip2);
        const _v4sf t1 = __builtin_ia32_unpckhps(ip0, ip2);
        const _v4sf t2 = __builtin_ia32_unpcklps(ip1, ip3);
        const _v4sf t3 = __builtin_ia32_unpckhps(ip1, ip3);

        const _v4sf ipx = __builtin_ia32_unpcklps(t0, t2);
        const _v4sf ipy = __builtin_ia32_unpckhps(t0, t2);
        const _v4sf ipz = __builtin_ia32_unpcklps(t1, t3);
        const _v4sf iph = __builtin_ia32_unpckhps(t1, t3);
        const _v4sf iph2 = iph*iph;

        _v4sf inb = {0.0f,0.0f,0.0f,0.0f};
        for (int j = 0; j < nj2; j += 2)
        {
          const _v4sf jp = *(jb + j);

#if 0 /*makes it slow*/
          const bool skip    = __builtin_ia32_movmskps(__builtin_ia32_orps(
                __builtin_ia32_cmpltps(jp,   imin),
                __builtin_ia32_cmpltps(imax, jp))) & 7;
          if (skip) continue;
#endif

          const _v4sf jpx = __builtin_ia32_shufps(jp, jp, 0x00);
          const _v4sf jpy = __builtin_ia32_shufps(jp, jp, 0x55);
          const _v4sf jpz = __builtin_ia32_shufps(jp, jp, 0xAA);


          const _v4sf dx = jpx - ipx;
          const _v4sf dy = jpy - ipy;
          const _v4sf dz = jpz - ipz;
          const _v4sf r2 = dx*dx + dy*dy + dz*dz;

          const _v4sf mask = __builtin_ia32_cmpltps(r2, iph2);
#if 0
          const int imask = __builtin_ia32_movmskps(mask);
          if (imask == 0) continue;
#endif
          inb += __builtin_ia32_andps((_v4sf){1,1,1,1}, mask);

        }
        nb[i+0] += (int)__builtin_ia32_vec_ext_v4sf(inb, 0);
        nb[i+1] += (int)__builtin_ia32_vec_ext_v4sf(inb, 1);
        nb[i+2] += (int)__builtin_ia32_vec_ext_v4sf(inb, 2);
        nb[i+3] += (int)__builtin_ia32_vec_ext_v4sf(inb, 3);
        asm("#eg02");
      }
#else
      for (int i = 0; i < igroup.nb(); i++)
      {
        const Body &ibody = igroup[i];
        const float4 ipos = ibody.pos_h();
        const real     h  = ibody.h();
        const real     h2 = h*h;
        if (overlapped(boundary(ipos), leafBnd))
          for (int j = 0; j < leaf.nb(); j++)
          {
            const float4 jpos = leaf[j].pos_h();
            const float  r2   = (ipos - jpos).norm2();
            if (r2 < h2)
              nb[i]++;
          }
      }
#endif
    }
  }
}

#endif /* __RANGE_SEARCH_GROUP_H__ */
