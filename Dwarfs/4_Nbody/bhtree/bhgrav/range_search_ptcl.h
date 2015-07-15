#ifndef _RANGE_SEARCH_PTCL_H__
#define _RANGE_SEARCH_PTCL_H__
  
template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
int range_search(const Body &ibody, const boundary &ibnd = boundary(), const int addr = 0, int nb = 0) const
{
  if (ROOT)
  {
    assert(isTreeReady());
    const boundary ibnd(ibody.pos_h());
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
        if (!not_overlapped(ibnd, bndsList[cellList[k].id()].inner()))
          nb = range_search<false>(ibody, ibnd, k, nb);
  }
  else
  {
    const Cell cell = cellList[addr];
    if (cell.isNode())
    {
      for (int k = 0; k < 8; k++)
        if (!cellList[cell.addr()+k].isEmpty())
          if (!not_overlapped(ibnd, bndsList[cellList[cell.addr()+k].id()].inner()))
            nb = range_search<false>(ibody, ibnd, cell.addr()+k, nb);
    }
    else
    {
      const Leaf &leaf  = leafList[cell.leafIdx()];
      const float4 ipos = ibody.pos_h();
      const real     h2 = ibody.h()*ibody.h();
      for (int i = 0; i < leaf.nb(); i++)
      {
        const float4 jpos = leaf[i].pos_h();
        const float  r2   = (ipos - jpos).norm2();
        if (r2 < h2)
          nb++;
      }
    }
  }

  return nb;
}

#endif /* _RANGE_SEARCH_PTCL_H__ */
