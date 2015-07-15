#ifndef __GROUP_H__
#define __GROUP_H__


template<const int N>
struct GroupT 
{
  typedef std::vector<GroupT> Vector;
  private:
  Body _list[N];     /* idx of the body */
  int _nb;           /* number of bodies in the list */
  int iPad[7];       /* alignment */

  public:
  GroupT() : _nb(0) {
    //        assert(8*(N+1)*sizeof(float) == sizeof(GroupT));
  }
  GroupT(const Body &body) : _nb(0) {insert(body);}
  void insert(const Body &body) 
  {
    assert(_nb < N);
    _list[_nb++] = body;
  }
  void remove(const int idx)
  {
    assert(idx < _nb);
    _nb--;
    for (int i = idx; i < _nb; i++)
      _list[i] = _list[i+1];
  }
  bool  isFull() const {return N == _nb;}
  bool isEmpty() const {return 0 == _nb;}
  const Body& operator[](const int i) const { return _list[i];}
  int nb() const {return _nb;}

  boundary innerBoundary() const
  {
    boundary bnd;
    for (int i = 0; i < _nb; i++)
      bnd.merge(_list[i].vector_pos());
    return bnd;
  }
  boundary outerBoundary() const
  {
    boundary bnd;
    for (int i = 0; i < _nb; i++)
      bnd.merge(_list[i].pos_h());
    return bnd;
  }
  Boundaries computeBoundaries() const
  {
    return Boundaries(innerBoundary(), outerBoundary());
  }

  boundary sort()   /* does peano-hilbert sort and returns outer boundary */
  {
    const boundary bnd = outerBoundary();
    const vec3    hlen = bnd.hlen();
    const float   size = __max(hlen.x, __max(hlen.y, hlen.z)) * 2.0f;
    const int     dfac = 1.0f / size * (PeanoHilbert::Key(1) << PeanoHilbert::BITS_PER_DIMENSION);
    const vec3     min = bnd.min;
    PeanoHilbert::PH  key_list[N];
    Body             body_list[N];
    for (int i = 0; i < _nb; i++)
    {
      key_list [i] = PeanoHilbert::PH(_list[i].vector_pos() - min, dfac, i);
      body_list[i] = _list[i];
    }
    std::sort(key_list, key_list+_nb, PeanoHilbert::PH());
    for (int i = 0; i < _nb; i++)
      _list[i] = body_list[key_list[i].idx];
    return bnd;
  }
} __attribute__ ((aligned(32)));

#endif /* __GROUP_H__ */
