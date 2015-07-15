#ifndef __MULTIPOLE_H__
#define __MULTIPOLE_H__

template<typename real>
struct Monopole
{
  typedef std::vector<Monopole> Vector;
  typedef vector3<real> vec3;
  private:
  vec3 _mpos;
  real _mass;

  public:

  Monopole() : _mpos(vec3((real)0.0)), _mass((real)0.0) {}
  Monopole(const vec3 &pos, const real mass) :
    _mpos(pos*mass), _mass(mass) {}

  template<typename fp>
    Monopole<real>(const Monopole<fp> &m) : _mpos(m.mpos()), _mass(m.mass()) {}

  const Monopole& operator+=(const Monopole &m)
  {
    _mpos += m._mpos;
    _mass += m._mass;
    return *this;
  }

  void CoM() 
  {
    assert(_mass > (real)0.0);
    _mpos = _mpos/_mass;
    _mass = _mass;
  }

  const vec3& mpos() const {return _mpos;}
  real        mass() const {return _mass;}
};

template<typename real>
struct Quadrupole
{
  typedef std::vector<Quadrupole> Vector;
  typedef vector3<real> vec3;
  enum TYPE {UNIT};

  private:
  real _xx, _yy, _zz, iPad1;
  real _xy, _xz, _yz, iPad2;

  public:
  Quadrupole() : _xx(0.0), _yy(0.0), _zz(0.0), _xy(0.0), _xz(0.0), _yz(0.0) {}
  Quadrupole(const real x, const TYPE type) : _xx(x), _yy(x), _zz(x), _xy(0.0), _xz(0.0), _yz(0.0) {}
  Quadrupole(const vec3 &pos, const real mass) :
    _xx(mass*pos.x*pos.x),
    _yy(mass*pos.y*pos.y),
    _zz(mass*pos.z*pos.z),
    _xy(mass*pos.x*pos.y),
    _xz(mass*pos.x*pos.z),
    _yz(mass*pos.y*pos.z) {}

  template<typename fp>
    Quadrupole<real>(const Quadrupole<fp> &q) : 
      _xx(q.xx()), _yy(q.yy()), _zz(q.zz()),
      _xy(q.xy()), _xz(q.xz()), _yz(q.yz()) {}

  const Quadrupole& operator+=(const Quadrupole &q)
  {
    _xx += q._xx; _yy += q._yy; _zz += q._zz;
    _xy += q._xy; _xz += q._xz; _yz += q._yz;
    return *this;
  }
  const Quadrupole& operator*=(const real x)
  {
    _xx *= x; _yy *= x; _zz *= x;
    _xy *= x; _xz *= x; _yz *= x;
    return *this;
  }

  real xx() const {return _xx;}
  real yy() const {return _yy;}
  real zz() const {return _zz;}
  real xy() const {return _xy;}
  real xz() const {return _xz;}
  real yz() const {return _yz;}
  real trace() const {return _xx + _yy + _zz;}

};

template<typename real>
struct Multipole
{
  typedef std::vector<Multipole> Vector;
  typedef vector3<real> vec3;
  typedef   Monopole<real> Mono;
  typedef Quadrupole<real> Quad;

  private:

  Mono   _monopole;  // +4= 4
  Quad _quadrupole;  // +8= 12

  public:

  Multipole(const Mono &m = Mono(), const Quad& q = Quad()) : _monopole(m), _quadrupole(q) 
  {
    assert(sizeof(Multipole) == sizeof(real)*12);
  }
  Multipole(const vec3 &pos, const real mass) : _monopole(pos, mass), _quadrupole(pos, mass) {}

  template<typename fp>
    Multipole<real>(const Multipole<fp> &m) : 
      _monopole(m.monopole()), _quadrupole(m.quadrupole()) {}

  const Multipole& operator+=(const Multipole &m)
  {
    _monopole   += m.  _monopole;
    _quadrupole += m._quadrupole;
    return *this;
  }

  const Mono   monopole() const {return   _monopole;}
  const Quad quadrupole() const {return _quadrupole;}

  Multipole complete() const 
  {
    Multipole m(*this);
    m._monopole.CoM();
    const real _mass = m.  _monopole. mass();
    const vec3  _pos = m.  _monopole. mpos();
    const real trace = m._quadrupole.trace();
    const real P     = trace - _mass*_pos.norm2();
    m._quadrupole += Quadrupole<real>(_pos, -_mass);
    m._quadrupole *= (real)3.0;
    m._quadrupole += Quadrupole<real>(-P, Quadrupole<real>::UNIT);


    return m;
  }
};

typedef Monopole<double> dMonopole;
typedef Monopole<float > fMonopole;

typedef Quadrupole<double> dQuadrupole;
typedef Quadrupole<float > fQuadrupole;

typedef Multipole<double> dMultipole;
typedef Multipole<float > fMultipole;

  template<const bool ROOT>
dMultipole computeMultipole(const Particle::Vector &ptcl, const int addr = 0)
{
  dMultipole multipole;
  if (ROOT)
  {
    assert(isTreeReady());
    multipoleList.resize(ncell);
    cellCoM      .resize(ncell);
    for (int k = 0; k < 8; k++)
      if (!cellList[k].isEmpty() && cellList[k].isTouched())
        multipole += computeMultipole<false>(ptcl, k);
    multipole = multipole.complete();
  }
  else
  {
    Cell &cell = cellList[addr];
    assert(!cell.isEmpty());
    if (cell.isNode())
    {
      const int i = cell.addr();
      for (int k = 0; k < 8; k++)
        if (!cellList[i+k].isEmpty() && cellList[i+k].isTouched())
          multipole += computeMultipole<false>(ptcl, i+k);
    }
    else
    {
      const Leaf &leaf = leafList[cell.leafIdx()];
      for (int i = 0; i < leaf.nb(); i++)
      {
        const int idx = leaf[i].idx();
        assert(idx < (int)ptcl.size());
        multipole += dMultipole(ptcl[idx].pos, ptcl[idx].mass);
        assert(overlapped(bndsList[cell.id()].inner(), leaf[i ].vector_pos()));
        assert(overlapped(bndsList[cell.id()].inner(), ptcl[idx].pos));
      }
    }
    multipoleList[cell.id()] = fMultipole(multipole.complete());

    /* compute opening criterion for this cell */

    const int id = cell.id();
    const fMultipole &m = multipoleList[id];
    const vec3 com = m.monopole().mpos();
    const float  s = (bndsList[id].inner().center() - com).abs();
    const vec3 len =  bndsList[id].inner().  hlen();
    const float  l = __max(len.x, __max(len.y, len.z)) * 2.0f;
    cellCoM[id] = float4(com.x, com.y, com.z, SQR(l*inv_theta + s));

#if 0
    fprintf(stderr, "dcom= %f %f %f  c= %f %f %f  l=%G\n",
        com.x, com.y, com.z,
        bndsList[id].inner().center().x,
        bndsList[id].inner().center().y,
        bndsList[id].inner().center().z, l);

    static int cnt = 0;
    if (cnt++ > 1000) exit(0);
#endif

    /* now the physical properties of the cell have been updated, unTouch it */

    cell.unsetTouched();
  }
  return multipole;
}

#endif /* __MULTIPOLE_H__ */
