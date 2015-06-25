#ifndef __VECTOR3_H
#define __VECTOR3_H

#include <iostream>
#include <fstream>
#include <cmath>

struct float4
{
#ifdef __SSE_H__
  _v4sf vec;
  float4(const _v4sf _vec) : vec(_vec){}
  float4(const float _x, const float _y, const float _z, const float _w) 
  {
    vec = (_v4sf){_x, _y, _z, _w};
  }
  float4(const float _x)
  {
    vec = (_v4sf){_x, _x, _x, _x};
  }
  operator _v4sf() const {return vec;}
  float4 operator-(const float4 rhs) const
  {
    return vec - (_v4sf)rhs;
  }
  float4 operator+(const float4 rhs) const
  {
    return (_v4sf)rhs + vec;
  }
  float4 operator*(const float4 rhs) const
  {
    return (_v4sf)rhs * vec;
  }
  float4& operator+=(const float4 v)
  {
    vec += v.vec;
    return *this;
  }
  float norm2() const
  {
    const _v4si mask = (_v4si){-1, -1, -1, 0}; 
    const _v4sf r    = __builtin_ia32_andps(vec, (_v4sf)mask);
    const _v4sf r2   = r*r;
    const _v4sf tmp  = __builtin_ia32_haddps(r2,  r2);
    const _v4sf res  = __builtin_ia32_haddps(tmp, tmp);
    return __builtin_ia32_vec_ext_v4sf(res, 0);
  }
  friend const float4 fabs(const float4 x)
  {
    const _v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
    return float4(__builtin_ia32_andps(x.vec, (_v4sf)mask));
  }
  float reduce() const
  {
    _v4sf v = __reduce_v4sf(vec);
    return __builtin_ia32_vec_ext_v4sf(v, 0);
  }

  /* accessor methods */

  template<int e>
    struct Accessor
    {
      _v4sf &vec;
      Accessor(_v4sf &_vec) : vec(_vec) {}
      operator float() 
      {
        return __builtin_ia32_vec_ext_v4sf(vec, e);
      }
      float4 operator=(const float f)
      {
        vec = __builtin_ia32_vec_set_v4sf(vec, f, e);
        return float4(vec);
      }
    };

  template<int e>
    struct ConstAccessor
    {
      const _v4sf &vec;
      ConstAccessor(const _v4sf &_vec) : vec(_vec) {}
      operator float() 
      {
        return __builtin_ia32_vec_ext_v4sf(vec, e);
      };
      private:
      float4 operator=(const float f) 
      {
        assert(0);
        return float4(vec);
      }
    };

  Accessor<0> x() {return Accessor<0>(vec);}
  Accessor<1> y() {return Accessor<1>(vec);}
  Accessor<2> z() {return Accessor<2>(vec);}
  Accessor<3> w() {return Accessor<3>(vec);}
  ConstAccessor<0> x() const {return ConstAccessor<0>(vec);}
  ConstAccessor<1> y() const {return ConstAccessor<1>(vec);}
  ConstAccessor<2> z() const {return ConstAccessor<2>(vec);}
  ConstAccessor<3> w() const {return ConstAccessor<3>(vec);}
#else
  float _x, _y, _z, _w;
  float4(const float x, const float y, const float z, const float w) 
  {
    _x = x;
    _y = y;
    _z = z;
    _w = w;
  }
  float4(const float x)
  {
    _x = x;
    _y = x;
    _z = x;
    _w = x;
  }
  float4 operator-(const float4 v) const
  {
    return float4(_x-v._x, _y-v._y, _z-v._z, _w-v._w);
  }
  float4 operator+(const float4 v) const
  {
    return float4(_x+v._x, _y+v._y, _z+v._z, _w+v._w);
  }
  float4 operator*(const float4 v) const
  {
    return float4(_x*v._x, _y*v._y, _z*v._z, _w*v._w);
  }
  float norm2() const
  {
    return _x*_x+_y*_y+_z*_z;
  }
  float4& operator+=(const float4 v)
  {
    _x += v._x;
    _y += v._y;
    _z += v._z;
    _w += v._w;
    return *this;
  }
  float& x() {return _x;}
  float& y() {return _y;}
  float& z() {return _z;}
  float& w() {return _w;}
  float x() const {return _x;}
  float y() const {return _y;}
  float z() const {return _z;}
  float w() const {return _w;}

  friend const float4 maxeach (const float4 &a, const float4 &b)
  {
    return float4(std::max(a._x, b._x), std::max(a._y, b._y), std::max(a._z, b._z), 0.0f);
  }
  friend const float4 mineach (const float4 &a, const float4 &b)
  {
    return float4(std::min(a._x, b._x), std::min(a._y, b._y), std::min(a._z, b._z), 0.0f);
  }
  friend const float4 fabs(const float4 &a)
  {
    return float4(std::abs(a._x), std::abs(a._y), std::abs(a._z), std::abs(a._w));
  }
  float reduce() const
  {
    return _x + _y + _z + _w;
  }
#endif
  float4() {}
};


template <class REAL> struct vector3{
  public:
    REAL x, y, z;
    vector3 (const float4 r) : x(r.x()), y(r.y()), z(r.z()) {}

    vector3(){
      x = y = z = REAL(0);
    }
    vector3(const REAL &r){
      x = y = z = r;
    }
    vector3(const REAL &_x, const REAL &_y, const REAL &_z){
      x = _x;  y = _y;  z = _z;
    }
    vector3(const REAL *p){
      x = p[0]; y = p[1]; z = p[2];
    }
    ~vector3(){}

    REAL &operator [](int i){
      return (&x)[i];
    }
    const REAL &operator [](int i) const{
      return (&x)[i];
    }
    template <class real> 
      operator vector3<real> () const {
        return vector3<real> (real(x), real(y), real(z));
      }
    operator REAL *(){
      return &x;
    }
    float4 to_float4() const {return float4(x,y,z,0.0f);}
    REAL (*toPointer())[3]{
      return (REAL (*)[3])&x;
    }
    typedef REAL (*pArrayOfReal3)[3];
    operator pArrayOfReal3(){
      return toPointer();
    }

    void outv(std::ostream &ofs = std::cout) const{
      ofs << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }
    bool are_numbers () const{
      // returns false if *this has (a) NaN member(s)
      return (norm2() >= REAL(0));
    }

    REAL norm2() const{
      return (*this)*(*this);
    }
    REAL abs() const{
      return std::sqrt(norm2());
    }

    friend std::ostream &operator << (std::ostream &ofs, const vector3<REAL> &v){
      ofs << v.x << " " << v.y << " " << v.z;
      return ofs;
    }
    friend std::istream &operator >> (std::istream &ifs, vector3<REAL> &v){
      ifs >> v.x >> v.y >> v.z;
      return ifs;
    }
    const vector3<REAL> operator + (const vector3<REAL> &v) const{
      return vector3<REAL> (x+v.x, y+v.y, z+v.z);
    }
    const inline vector3<REAL> operator - (const vector3<REAL> &v) const{
      return vector3<REAL> (x-v.x, y-v.y, z-v.z);
    }
    const vector3<REAL> operator * (const REAL &s) const{
      return vector3<REAL> (x*s, y*s, z*s);
    }
    friend const vector3<REAL> operator * (const REAL &s, const vector3<REAL> &v){
      return v*s;
    }
    // dot product
    const inline REAL operator * (const vector3<REAL> &v) const{
      return (x*v.x + y*v.y + z*v.z);
    }
    // vector product
    const vector3<REAL> operator % (const vector3<REAL> &v) const{
      return vector3<REAL> (
          y*v.z - z*v.y, 
          z*v.x - x*v.z, 
          x*v.y - y*v.x);
    }
    const vector3<REAL> operator / (const REAL &s) const{
      REAL r = REAL(1)/s;
      return (*this)*r;
    }
    const vector3<REAL> operator = (const vector3<REAL> &v){
      x = v.x; y=v.y; z=v.z;
      return *this;
    }

    const vector3<REAL> operator - (){
      return vector3<REAL> (-x, -y, -z);
    }
    const vector3<REAL> &operator += (const vector3<REAL> &v){
      *this = *this + v;
      return *this;
    }
    const vector3<REAL> &operator -= (const vector3<REAL> &v){
      *this = *this - v;
      return *this;
    }
    const vector3<REAL> &operator *= (const REAL &s){
      *this = *this * s;
      return *this;
    }
    const vector3<REAL> &operator /= (const REAL &s){
      *this = *this / s;
      return *this;
    }

    friend const vector3<REAL> maxeach (const vector3<REAL> &a, const vector3<REAL> &b){
      return vector3<REAL> (std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
    }
    friend const vector3<REAL> mineach (const vector3<REAL> &a, const vector3<REAL> &b){
      return vector3<REAL> (std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
    }
    const vector3<REAL> abseach(){
      return vector3<REAL> (std::fabs(x), std::fabs(y), std::fabs(z));
    }
};

typedef vector3<double> dvec3;
typedef vector3<float>  fvec3;



#endif //  __VECTOR3_H
