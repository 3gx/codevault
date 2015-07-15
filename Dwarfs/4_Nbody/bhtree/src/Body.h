#ifndef __BODY_H__
#define __BODY_H__
  
struct Body
{
  typedef std::vector<Body> Vector;
  private:
  float4 packed_pos;  /* x y z h */
  float4 packed_mass;  /* idx iPad1 iPad2 mass */

  public:
  Body() 
  {
    assert(sizeof(Body) == sizeof(float)*8);
  }
  Body(const vec3 &pos, const int idx, const float h = 0.0f, const float mass = 0.0f) :
    packed_pos(pos.x, pos.y, pos.z, h), packed_mass((float)idx, 0.0f, 0.0f, mass) {}
  Body(const Particle &p, const int idx) : 
    packed_pos(p.pos.x, p.pos.y, p.pos.z, p.h), packed_mass((float)idx, 0.0f, 0.0f, p.mass) {}
  Body(const vec3 &pos, const float h, const int idx = -1, const float mass = 0.0f) : 
    packed_pos(pos.x, pos.y, pos.z, h), packed_mass((float)idx, 0.0f, 0.0f, mass) {}
  int           idx() const {return (int)packed_mass.x();}
  float        mass() const {return      packed_mass.w();}
  float           h() const {return      packed_pos .w();}
  float4      pos_h() const {return      packed_pos;     }
  vec3   vector_pos() const {return vec3(packed_pos.x(), packed_pos.y(), packed_pos.z());}
#ifdef __SSE_H__
  Body& operator=(const Body &rhs)
  {
    _v4sf *dst = (_v4sf*)this;
    _v4sf *src = (_v4sf*)&rhs;
#if 1
    dst[0] = src[0];
    dst[1] = src[1];
#else
    const v4sf a = src[0];
    const v4sf b = src[1];
    __builtin_ia32_movntps((float*)(dst+0), a);
    __builtin_ia32_movntps((float*)(dst+1), b);
#endif
    return *this;
  }

  float4 pos_mass() const 
  {
    const _v4sf res = __builtin_ia32_blendps(packed_pos, packed_mass, 1<<3);
    return  res;
  }
#else
  float4   pos_mass() const 
  {
    return float4(packed_pos.x(), packed_pos.y(), packed_pos.z(), packed_mass.w()); 
  }
#endif
};

#endif /* __BODY_H__ */
