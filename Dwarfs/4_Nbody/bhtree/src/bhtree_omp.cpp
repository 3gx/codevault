#include <cstdio>
#include <cstdlib>
#include "octree.h"
#include "plummer.h"
#include "mytimer.h"

enum {NGROUP = 96};
typedef Octree::GroupT<NGROUP> octGroup;

int main(int argc, char * argv[])
{
  int n_bodies = 10240;
  real theta = 0.75;
  if (argc > 1) n_bodies = atoi(argv[1]);
  if (argc > 2) theta = atof(argv[2]);
  assert(n_bodies > 0);
  fprintf(stderr, "n_bodies= %d\n", n_bodies);
  fprintf(stderr, "theta= %g\n", theta);

  const double t00 = get_wtime();
  Particle::Vector ptcl;
  ptcl.reserve(n_bodies);
#if 1
#define PLUMMER
  const Plummer data(n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
    ptcl.push_back(Particle(data.mass[i], data.pos[i], data.vel[i], i));
  }
#elif 0  /* reads IC */
  int dummy;
  std::cin >> dummy >> n_bodies;
  fprintf(stderr, " -- input file: nbodies= %d\n", n_bodies);
  for (int i = 0; i < n_bodies; i++)
  {
    Particle p;
    std::cin >> p.pos.x >> p.pos.y >> p.pos.z >> p.h >> p.nb;
    p.h *= 2.0;
    ptcl.push_back(p);
  }
#elif 1
  {
    const int nb_mean = 128;
    const real s = std::pow(3.0/(4.0*M_PI)*(double)nb_mean/(double)n_bodies, 1.0/3.0);
    for (int i = 0; i < n_bodies; i++)
    {
      ptcl.push_back(Particle(1.0/n_bodies,
            vec3(drand48(), drand48(), drand48()),
            vec3(0.0, 0.0, 0.0),
            1, s));
    }
  }
#endif

  Octree::Body::Vector octBodies;
  octBodies.reserve(n_bodies);

  const double t10 = get_wtime();
  fprintf(stderr, " -- Compute bounding box -- \n");

  vec3 rmin(+HUGE);
  vec3 rmax(-HUGE);

  for (int i = 0; i < n_bodies; i++)
  {
    octBodies.push_back(Octree::Body(ptcl[i], i));
    rmin = mineach(rmin, ptcl[i].pos);
    rmax = maxeach(rmax, ptcl[i].pos);
  }
  const vec3 centre = (rmax + rmin)*0.5;
  const vec3 vsize  =  rmax - rmin;
  const real  size  = __max(__max(vsize.x, vsize.y), vsize.z);
  real size2 = 1.0;
  while (size2 > size) size2 *= 0.5;
  while (size2 < size) size2 *= 2.0;

  const int n_nodes = n_bodies;
  Octree tree(centre, size2, n_nodes, theta);

  const double t20 = get_wtime();
  fprintf(stderr, " -- Buidling octTree -- \n");
  for (int i = 0; i < n_bodies; i++)
  {
    tree.insert(octBodies[i]);
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d np= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth(), tree.get_np());
  const double t30 = get_wtime();

  fprintf(stderr, " -- Dump morton -- \n");
  std::vector<int> morton_list;
  morton_list.reserve(n_bodies);
  tree.tree_dump<true>(morton_list);
  fprintf(stderr, "morton_list.size()= %d  n_bodies= %d\n",
      (int)morton_list.size(), n_bodies);
  for (std::vector<int>::iterator it = morton_list.begin(); it != morton_list.end(); it++)
    assert(*it < n_bodies);
  assert((int)morton_list.size() == n_bodies);

  const double t40 = get_wtime();
  fprintf(stderr, " -- Shuffle octBodies -- \n");
  Particle::Vector sortedPtcl;
  sortedPtcl.reserve(n_bodies);
  for (std::vector<int>::iterator it = morton_list.begin(); it != morton_list.end(); it++)
  {
    assert(*it < n_bodies);
    sortedPtcl.push_back(ptcl[octBodies[*it].idx()]);
  }

  const double t50 = get_wtime();
  fprintf(stderr, " -- Buidling octTreeSorted -- \n");
  tree.clear();
  for (int i = 0; i < n_bodies; i++)
  {
    tree.insert(Octree::Body(sortedPtcl[i], i));
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());

  const double t60 = get_wtime();
  fprintf(stderr, " -- Compute boundaries -- \n");
  tree.computeBoundaries<true>();
  const boundary root_innerBnd = tree.root_innerBoundary();
  fprintf(stderr, " rootBnd_inner= %g %g %g  size= %g %g %g \n",
      root_innerBnd.center().x,
      root_innerBnd.center().y,
      root_innerBnd.center().z,
      root_innerBnd.hlen().x,
      root_innerBnd.hlen().y,
      root_innerBnd.hlen().z);
  const boundary root_outerBnd = tree.root_outerBoundary();
  fprintf(stderr, " rootBnd_outer= %g %g %g  size= %g %g %g \n",
      root_outerBnd.center().x,
      root_outerBnd.center().y,
      root_outerBnd.center().z,
      root_outerBnd.hlen().x,
      root_outerBnd.hlen().y,
      root_outerBnd.hlen().z);
  fprintf(stderr, "rootCentre:= %g %g %g  rootSize= %g \n",
      tree.get_rootCentre().x,
      tree.get_rootCentre().y,
      tree.get_rootCentre().z,
      tree.get_rootSize());

  const double t63 = get_wtime();
  assert(tree.sanity_check<true>() == n_bodies);
  const double t65 = get_wtime();
  tree.buildLeafList<true>();
  const double t68 = get_wtime();

  octGroup::Vector groupList;
  groupList.reserve(tree.nLeaf());

  const bool SORT = 0 ? true : false;  /* use peano-sort inside the group */
  tree.buildGroupList<SORT, true>(groupList);

  const double t69 = get_wtime();

  fprintf(stderr, " -- Range search -- \n");
  int nb = 0;
#ifndef PLUMMER
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
#if 0
    nb += tree.range_search<true>(octBodies[i]);
#else
    nb += tree.range_search<true>(Octree::Body(sortedPtcl[i], i));
#endif
  }
#endif
  const double t70 = get_wtime();
  const int ngroup = groupList.size();
  fprintf(stderr, " -- Range search Group-Leaf : ngroup=%d  nbody= %d-- \n", ngroup, n_bodies);
  int nbL = 0;
#ifndef PLUMMER
#pragma omp parallel for reduction(+:nbL)
  for (int i = 0; i < ngroup; i++)
  {
    int nb[NGROUP];
    const octGroup &group = groupList[i];
    tree.range_search<true>(nb, group);
    for (int j = 0; j < group.nb(); j++)
    {
#if 1
      const int idx = group[j].idx();
      const int nb0 = ptcl[idx].nb;
      if (nb0 != nb[j] && nb0 > 0)
        fprintf(stderr, "nb0= %d  nb= %d\n", nb0, nb[j]);
#endif
      nbL += nb[j];
    }
  }
#endif
  const double t75 = get_wtime();

#if 0
  fprintf(stderr, " -- Remove ptcl -- \n");
  int nrm = 0;
#if 1
  for (int i = 0; i < n_bodies; i++)
  {
    tree.remove(Octree::Body(ptcl[i].pos, i));
    nrm++;
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());
#endif

  const double t80 = get_wtime();

  fprintf(stderr, " -- Insert ptcl -- \n");
  int nins = 0;
#if 1
  for (int i = 0; i < n_bodies; i++)
  {
    tree.insert(octBodiesSorted[i]);
    nins++;
  }
  fprintf(stderr, "ncell= %d nnode= %d nleaf= %d n_nodes= %d  depth= %d\n",
      tree.get_ncell(), tree.get_nnode(), tree.get_nleaf(), n_nodes, tree.get_depth());
#endif

  const double t90 = get_wtime();

  fprintf(stderr, " -- Range search 1 -- \n");
  int nb1 = 0;
#if 1
#pragma omp parallel for reduction(+:nb)
  for (int i = 0; i < n_bodies; i++)
  {
#if 0
    nb1 += tree.range_search<true>(octBodies[i]);
#else
    nb1 += tree.range_search<true>(octBodiesSorted[i]);
#endif
  }
#endif
#endif
  const double t100 = get_wtime();

#if 1
  {
    fprintf(stderr, " Computing multipole \n");
    const Octree::dMultipole rootM = tree.computeMultipole<true>(sortedPtcl);
    fprintf(stderr,  " Mass= %g \n", rootM.monopole().mass());
    const vec3 mpos = rootM.monopole().mpos();
    fprintf(stderr, " Monopole= %g %g %g  \n", mpos.x, mpos.y, mpos.z);
    fprintf(stderr, " Quadrupole: xx= %g yy= %g zz= %g trace= %g  xy= %g xz= %g zz= %g \n",
        rootM.quadrupole().xx(),
        rootM.quadrupole().yy(),
        rootM.quadrupole().zz(),
        rootM.quadrupole().trace(),
        rootM.quadrupole().xy(),
        rootM.quadrupole().xz(),
        rootM.quadrupole().yz());
  }
#endif
  double t110 = get_wtime();
  {
    Octree::dMultipole M;
    for (int i = 0; i < n_bodies; i++)
      M += Octree::dMultipole(sortedPtcl[i].pos, ptcl[i].mass);
    M = M.complete();
    fprintf(stderr,  " Mass= %g \n", M.monopole().mass());
    const vec3 mpos = M.monopole().mpos();
    fprintf(stderr, " Monopole= %g %g %g  \n", mpos.x, mpos.y, mpos.z);
    fprintf(stderr, " Quadrupole: xx= %g yy= %g zz= %g trace= %g  xy= %g xz= %g zz= %g \n",
        M.quadrupole().xx(),
        M.quadrupole().yy(),
        M.quadrupole().zz(),
        M.quadrupole().trace(),
        M.quadrupole().xy(),
        M.quadrupole().xz(),
        M.quadrupole().yz());
  }
  t110 = get_wtime();

#if 1
  {
    fprintf(stderr, " Computing gForce \n");
    real   gpot = 0.0;
    double mtot = 0.0;
    double fx = 0, fy = 0, fz = 0;
    unsigned long long npc = 0, npp = 0;
#pragma omp parallel for reduction(+:mtot, gpot, fx, fy, fz, npc, npp)
    for (int i = 0; i < ngroup; i++)
    {
      const octGroup &group = groupList[i];
      float4 force[NGROUP] __attribute__ ((aligned(64)));

      const std::pair<int, int> ninter = tree.gForce(group, force);
      npp += ninter.first;
      npc += ninter.second;

      for (int j = 0; j < group.nb(); j++)
      {
        const int idx = group[j].idx();
        const Particle &p = ptcl[idx];
        mtot += p.mass;
        fx += p.mass * force[j].x();
        fy += p.mass * force[j].y();
        fz += p.mass * force[j].z();
        gpot += 0.5*p.mass*force[j].w();
      }
    }
    assert(mtot > 0.0);
    fprintf(stderr , "<Ng>= %g\n", (double)n_bodies/ngroup);
    fprintf(stderr, " Npp= %g  Npc= %g\n",
        (double)npp/n_bodies, (double)npc/n_bodies);
    fprintf(stderr, " mtot= %g  ftot= %g %g %g  gpot= %g\n",
        mtot, fx/mtot, fy/mtot, fz/mtot, gpot/mtot);
  }
#endif
  double t120 = get_wtime();


  fprintf(stderr, " Timing info: \n");
  fprintf(stderr, " -------------\n");
  fprintf(stderr, "   Plummer:  %g sec \n", t10 -t00);
  fprintf(stderr, "   BBox:     %g sec \n", t20 -t10);
  fprintf(stderr, "   Tree:     %g sec \n", t30 -t20);
  fprintf(stderr, "   Morton:   %g sec \n", t40 -t30);
  fprintf(stderr, "   Shuffle:  %g sec \n", t50 -t40);
  fprintf(stderr, "   TreeSort: %g sec \n", t60 -t50);
  fprintf(stderr, "   Boundary: %g sec \n", t63 -t60);
  fprintf(stderr, "   Sanity:   %g sec \n", t65 -t63);
  fprintf(stderr, "   LeafList: %g sec \n", t68 -t65);
  fprintf(stderr, "   GroupLst: %g sec \n", t69 -t68);
  fprintf(stderr, "   RangeS:   %g sec <nb>= %g \n", t70 -t69, (real)nb /n_bodies);
  fprintf(stderr, "   RangeL:   %g sec <nb>= %g \n", t75 -t70, (real)nbL/n_bodies);
#if 0
  fprintf(stderr, "   Remove:   %g sec nrm= %d \n", t80 - t70, nrm);
  fprintf(stderr, "   Insert:   %g sec nins= %d \n", t90 - t80, nins);
  fprintf(stderr, "   RangeS1:   %g sec <nb>= %g \n", t100 -t90, (real)nb1/n_bodies);
#endif
  fprintf(stderr, "   MultiP:   %g sec \n", t110 -t100);
  fprintf(stderr, "   gForce:   %g sec \n", t120 -t110);

  return 0;
}
