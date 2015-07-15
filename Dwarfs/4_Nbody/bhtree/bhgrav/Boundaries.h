#ifndef __BOUNDARIES_H__
#define __BOUNDARIES_H__

struct Boundaries
{
  typedef std::vector<Boundaries> Vector;
  private:
  boundary _inner;
  boundary _outer;

  public:
  Boundaries(const boundary &inner = boundary(), const boundary &outer = boundary()) :
    _inner(inner), _outer(outer) {}
  Boundaries& merge(const Boundaries &rhs) 
  {
    _inner.merge(rhs._inner);
    _outer.merge(rhs._outer);
    return *this;
  }

  const boundary& inner() const { return _inner; }
  const boundary& outer() const { return _outer; }

};

#endif /* __BOUNDARIES_H__ */
