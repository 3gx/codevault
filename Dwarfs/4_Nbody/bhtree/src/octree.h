#ifndef __OCTREE_H__
#define __OCTREE_H__

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <stack>
#include <algorithm>

#ifdef __mySSE__
#include "sse.h"
#endif

#ifdef __myAVX__
#include "avx.h"
#endif

#include "primitives.h"
#include "peano.h"
#include "boundary4.h"


struct Particle
{
  typedef std::vector<Particle> Vector;
  real mass;
  vec3 pos;
  vec3 vel;
  int  id;
  real h;
  int nb;

  Particle() {}
  Particle(const real _mass, const vec3 &_pos, const vec3 &_vel, const int _id, const real _h = 0.0) :
    mass(_mass), pos(_pos), vel(_vel), id(_id), h(_h), nb(0) {}
};


struct Octree
{
  enum {NLEAF =  16};

  private:
  enum {EMPTY =  -1};
  enum {BODYX =  -2};
  
  static inline int reverse_int(const int a) {return BODYX-a;}

  public:

#include "Body.h"
#include "Cell.h"
#include "Boundaries.h"
#include "Group.h"
#include "Multipole.h"

  typedef GroupT<NLEAF> Leaf;

  private:

  float4          root_centre;    /* root's geometric centre & size */
  int                   depth;    /* tree-depth */
  int                   nnode;    /* number of tree nodes */
  int                   ncell;    /* number of tree cells: node + leaf */
  bool              treeReady;    /* this is a flag if tree is ready for walks */
  Cell    ::Vector   cellList;    /* list of tree cells (node+leaf) */
  Leaf    ::Vector   leafList;    /* list of leaves */
  std::stack<int>    leafPool;    /* pool of available leaves */
  std::stack<int>    cellPool;    /* pool of available cells */
  Boundaries::Vector bndsList;    /* list of cell  boundaries */

  real inv_theta;                 /* 1/theta opening criterion */
  real eps2;                      /* softening squared */
  std::vector<float4> cellCoM; /* precomputed criterio tree-cells */

  fMultipole::Vector multipoleList;

  std::vector<int> leafList_addr;

  public:

  Octree(const vec3 &_centre, const real _size, const int n_nodes, const real _theta = 0.75, const real eps = 0.1) :
    depth(0), nnode(0), ncell(0), treeReady(false), inv_theta(1.0/__max(_theta, (real)1.0e-15)), eps2(eps*eps)
  {
    root_centre = float4(_centre.x, _centre.y, _centre.z, _size);
    cellList.resize(n_nodes<<3);
  }

  int get_depth() const { return depth;}
  int get_nnode() const { return nnode;}
  int get_ncell() const { return ncell;}
  int get_nleaf() const { return leafList.size();}

  vec3 get_rootCentre() const { return vec3(root_centre.x(), root_centre.y(), root_centre.z());}
  real get_rootSize  () const { return root_centre.w();}

  void clear(const int n_nodes = 0)
  {
    depth = nnode = ncell = 0;
    treeReady = false;
    const int nsize = n_nodes <= 0 ? cellList.size() : n_nodes << 3;
    cellList.clear();
    cellList.resize(nsize);
    leafList.clear();
    bndsList.clear();
    leafList_addr.clear();
    cellCoM.clear();
    leafPool = std::stack<int>();
    multipoleList.clear();
    cellCoM.clear();
  }
  bool isTreeReady() const {return treeReady;}


  inline int Octant(const float4 lhs, const float4 rhs) 
  {
#ifdef __SSE_H__
    int mask = __builtin_ia32_movmskps(
        __builtin_ia32_cmpgeps(rhs, lhs));
    return 7 & mask;
#else
    return
      (((lhs.x() <= rhs.x()) ? 1 : 0) +  
       ((lhs.y() <= rhs.y()) ? 2 : 0) + 
       ((lhs.z() <= rhs.z()) ? 4 : 0));
#endif
  }

  inline float4 child_centre(const float4 centre, const float4 ppos)
  {
#ifdef __SSE_H__
    const _v4si mask = {(int)0x80000000, (int)0x80000000, (int)0x80000000, 0};
    const _v4sf off  = {0.25f, 0.25f, 0.25f, 0.5f};
    const _v4sf len  = __builtin_ia32_shufps(centre, centre, 0xff);

    _v4sf tmp = __builtin_ia32_cmpgeps(ppos, centre); // mask bits
    tmp = __builtin_ia32_andps(tmp, (_v4sf)mask);    // sign mask
    tmp = __builtin_ia32_orps(tmp, off*len);        // offset;
    return (_v4sf)centre - tmp;
#else
    const int oct = Octant(centre, ppos);
    const float s = centre.w()*0.25f;
    return float4(
        centre.x() + s * ((oct&1) ? (real)1.0 : (real)-1.0),
        centre.y() + s * ((oct&2) ? (real)1.0 : (real)-1.0),
        centre.z() + s * ((oct&4) ? (real)1.0 : (real)-1.0),
        centre.w()*0.5f
        );
#endif
  }


  /**********************************/
  /*** Tree management primitives ***/
  /**********************************/

  void deleteCell(const int cellId)
  {
    cellPool.push(cellId);
  }
  void deleteLeaf(const int leafIdx)
  {
    leafPool.push(leafIdx);
  }
  template<const bool WITHBODY>
    Cell newLeaf(const Body &body = Body())
    {
      const int cell = cellPool.empty() ? ncell++ : cellPool.top();
      if (!cellPool.empty()) cellPool.pop();

      if (leafPool.empty())
      {
        const int addr = leafList.size();
        leafList.push_back(WITHBODY ? Leaf(body) : Leaf());
        return Cell(reverse_int(addr), cell, WITHBODY ? 1 : 0);
      }
      else
      {
        const int addr = leafPool.top();
        leafPool.pop();
        leafList[addr] = WITHBODY ? Leaf(body) : Leaf();
        return Cell(reverse_int(addr), cell, WITHBODY ? 1 : 0);
      }
    }

  int get_np() const 
  {
    int np = 0;
    for (int k = 0; k < 8; k++)
      np += cellList[k].np();
    return np;
  }
  
  /********************************************/
  /*** Tree management (insert/remove body) ***/
  /********************************************/

  void insert(const Body &new_body)
  {
    int  child_idx = 0;                      /* child idx inside a node */
    Cell child     = Cell(Cell::ROOT);       /* child */
    int locked     = 0;                      /* cell that needs to be updated */
    int depth      = 0;                      /* depth */ 
#ifndef NDEBUG
    const int n_nodes_max = cellList.size();
#endif

    float4 centre =  root_centre;

    /* walk the tree to find the first leaf or empty cell */
    while (child.isNode())  /* if non-negative, means it is a tree-node */
    {
      const Cell node = child;
      depth++;

      child_idx  =       Octant(centre, new_body.pos_h());
      centre     = child_centre(centre, new_body.pos_h());

      locked = node.addr() + child_idx;
      assert(locked < n_nodes_max);
      child  = cellList[locked];
      cellList[locked]++;
    }

    /* locked on the cell that needs to be updated */

    if (child.isEmpty())  /* the cell is empty, make it a leaf */
    {
      assert(child.isClean());
      cellList[locked] = newLeaf<true>(new_body);
      assert(cellList[locked].np() == leafList[cellList[locked].leafIdx()].nb());
    }
    else          /* this is already a leaf */
    {
      assert(child.isLeaf());
      int leaf_idx = child.leafIdx();
      assert(leaf_idx < (int)leafList.size());
      while(leafList[leaf_idx].isFull())   /* if leaf is full split it */
      {
        depth++;
        nnode++;
        const int cfirst = nnode<<3;
        assert(cfirst+7 < n_nodes_max);
        cellList[locked].set_addr(cfirst);

        assert(leaf_idx >= 0);
        assert(leaf_idx < get_nleaf());
        const Leaf leaf = leafList[leaf_idx];
        deleteLeaf(leaf_idx);

        for (int i = 0; i < leaf.nb(); i++)
        {
          const Body &body = leaf[i];
          child_idx  = Octant(centre, body.pos_h());
          Cell &cell = cellList[cfirst + child_idx];
          cell++;
          if (cell.isEmpty())
          {
            cell = newLeaf<true>(body);
          }
          else
          {
            assert(cell.isLeaf());
            assert(!leafList[cell.leafIdx()].isFull());
            leafList[cell.leafIdx()].insert(body);
          }
        }

        child_idx =       Octant(centre, new_body.pos_h());
        centre    = child_centre(centre, new_body.pos_h());

        locked = cfirst + child_idx;
        Cell &child  = cellList[locked];
        if (child.isEmpty())
          child = newLeaf<false>();
        child++;
        assert(child.isLeaf());
        leaf_idx = child.leafIdx();
      }
      leafList[leaf_idx].insert(new_body);  /* push particle into a leaf */
      assert(leafList[leaf_idx].nb() == cellList[locked].np());
    }

    this->depth = __max(this->depth, depth);
  }

  /**************/

  bool remove(const Body &body)
  {
    float4 centre =  root_centre;

#ifndef NDEBUG
    const int n_nodes_max = cellList.size();
#endif
    std::stack<int> path;

    Cell child     = Cell(Cell::ROOT);       /* child */
    int locked     = 0;                      /* cell that needs to be updated */
    int child_idx  = 0;
    while (child.isNode()) 
    {
      const Cell node = child;
      depth++;

      child_idx  =       Octant(centre, body.pos_h());
      centre     = child_centre(centre, body.pos_h());

      locked = node.addr() + child_idx;
      assert(locked < n_nodes_max);
      child  = cellList[locked];
      path.push(locked);
      cellList[locked].setTouched();
    }

    assert(!child.isEmpty());
    if (child.isEmpty()) return false;

    assert(child.isLeaf());
    const int leaf_idx = child.leafIdx();
    assert(leaf_idx < (int)leafList.size());
    Leaf &leaf = leafList[leaf_idx];
    int idx = 0;
    for (idx = 0; idx < leaf.nb(); idx++)
      if (body.idx() == leaf[idx].idx())
        break;

    assert(idx < leaf.nb());
    if (idx == leaf.nb()) return false;
    leaf.remove(idx);

    if (leaf.isEmpty()) 
    {
      while(!path.empty())
      {
        path.pop();
      }
      deleteLeaf(leaf_idx);
    }

    return true;
  } 

  /**************/

  template<const bool ROOT>  /* must be ROOT = true on the root node (first call) */
    void tree_dump(std::vector<int> &list, const int node = 0) const
    {
      if (ROOT)
      {
        for (int k = 0; k < 8; k++)
          tree_dump<false>(list, k);
      }
      else
      {
        assert(node < (int)cellList.size());
        const Cell cell = cellList[node];
        if (cell.isEmpty()) return;
        if (cell.isNode())
        {
          for (int k = 0; k < 8; k++)
            tree_dump<false>(list, cell.addr()+k);
        }
        else
        {
          assert(cell.isLeaf());
          const int leaf_idx = cell.leafIdx();
          assert(leaf_idx < get_nleaf());
          const Leaf &leaf = leafList[leaf_idx];
          for (int i = 0; i < leaf.nb(); i++)
            list.push_back(leaf[i].idx());
        }
      }
    }

  /**************/

  boundary root_innerBoundary()
  {
    boundary bnd;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
      {
        bnd.merge(bndsList[cellList[k].id()].inner());
      }
    return bnd;
  }
  boundary root_outerBoundary()
  {
    boundary bnd;
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty())
      {
        bnd.merge(bndsList[cellList[k].id()].outer());
      }
    return bnd;
  }

  template<const bool ROOT=true>  /* must be ROOT = true on the root node (first call) */
    Boundaries computeBoundaries(const int addr = 0)
    {
      Boundaries bnds;
      if (ROOT)
      {
        bndsList.resize(ncell);
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty() && cellList[k].isTouched())
          {
            computeBoundaries<false>(k);
            bnds.merge(bndsList[cellList[k].id()]);
          }
        treeReady = true;
      }
      else
      {
        assert(addr < (int)cellList.size());
        Cell &cell = cellList[addr];
        assert (!cell.isEmpty());
        assert(cell.id() >= 0);


        if (cell.isNode())
        {
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            if (!cellList[k].isEmpty() && cellList[k].isTouched())
            {
              computeBoundaries<false>(k);
              bnds.merge(bndsList[cellList[k].id()]);
            }
        }
        else
        {
          const Leaf &leaf = leafList[cell.leafIdx()];
          bnds.merge(leaf.computeBoundaries());
        }
        bndsList[cell.id()] = bnds;
//        cell.unsetTouched();
      }
      return bnds;
    }

  /**************/

  template<const bool ROOT=true>  /* must be ROOT = true on the root node (first call) */
    int sanity_check(const int addr = 0) const
    {
      int nb = 0;
      if (ROOT)
      {
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            nb += sanity_check<false>(k);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return nb;
        assert(cell.id() >= 0);

        if (cell.isNode())
        {
          assert(cell.np() > NLEAF);
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            nb += sanity_check<false>(k);
        }
        else
        {
          assert(cell.np() <= NLEAF);
          assert(cell.isLeaf());
          assert(cell.leafIdx() < get_nleaf());
          const Leaf &leaf = leafList[cell.leafIdx()];
          for (int i = 0; i < leaf.nb(); i++)
          {
            const vec3 jpos = leaf[i].vector_pos();
            assert(overlapped(bndsList[cell.id()].inner(), jpos));
            assert(overlapped(bndsList[cell.id()].outer(), leaf[i].pos_h()));
            nb++;
          }
        }
      }
      return nb;
    }

  /**************/

  int nLeaf() const {return leafList_addr.size();}
  const Leaf& getLeaf(const int i) const {return leafList[i];}
  template<const bool ROOT=true>  /* must be ROOT = true on the root node (first call) */
    void buildLeafList(const int addr = 0) 
    {
      if (ROOT)
      {
        leafList_addr.clear();
        leafList_addr.reserve(128);
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            buildLeafList<false>(k);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return;
        assert(cell.id() >= 0);

        if (cell.isNode())
        {
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            buildLeafList<false>(k);
        }
        else
        {
          assert(cell.isLeaf());
          assert(cell.leafIdx() < get_nleaf());
          leafList_addr.push_back(cell.leafIdx());
        }
      }
    }

  /**************/

  template<const bool SORT, const bool ROOT=true, const int N>  /* must be ROOT = true on the root node (first call) */
    void buildGroupList(std::vector< GroupT<N> > &groupList, 
        const int addr = 0) const
    {
      assert(NLEAF <= N);
      if (ROOT)
      {
        assert(isTreeReady());
        for (int k = 0; k < 8; k++)
          if (!cellList[k].isEmpty())
            buildGroupList<SORT, false>(groupList, k);
      }
      else
      {
        assert(addr < (int)cellList.size());
        const Cell &cell = cellList[addr];
        if (cell.isEmpty()) return;
        assert(cell.id() >= 0);

        if (cell.np() <= N)
        {
          groupList.push_back(GroupT<N>());
          extractBodies(groupList.back(), addr);
          if (SORT)
            groupList.back().sort();
        }
        else
        {
          assert(!cell.isLeaf());
          for (int k = cell.addr(); k < cell.addr()+8; k++)
            buildGroupList<SORT, false>(groupList, k);
        }
      }
    }

  template<const int N>
    void extractBodies(GroupT<N> &group, const int addr) const
    {
      assert(addr < (int)cellList.size());
      const Cell &cell = cellList[addr];
      if (cell.isEmpty()) return;

      if (cell.isNode())
      {
        for (int k = 0; k < 8; k++)
          extractBodies(group, cell.addr() + k);
      }
      else
      {
        assert(cell.isLeaf());
        const Leaf &leaf = leafList[cell.leafIdx()];
        for (int i = 0; i < leaf.nb(); i++)
          group.insert(leaf[i]);
      }
    }

  /**************/

#include "range_search_ptcl.h"
#include "range_search_group.h"
  
  /**************/

#include "gforce_group.h"

};


#endif /* __OCTREE_H__ */
