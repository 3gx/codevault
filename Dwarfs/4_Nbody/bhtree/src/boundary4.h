#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include <vector>
#include "vector3.h"
#include <cmath>

struct boundary{
  typedef std::vector<boundary> Vector;
	typedef vector3<float> vec3;

  float4 min, max;

	boundary() : min(HUGE), max(-HUGE) {}
  boundary(const float4 pos) 
  {
    min = float4(pos.x() - pos.w(), pos.y() - pos.w(), pos.z() - pos.w(), 0.0f);
    max = float4(pos.x() + pos.w(), pos.y() + pos.w(), pos.z() + pos.w(), 0.0f);
  }
	boundary(const float4 &_min, const float4 &_max) : min(_min), max(_max) {}
	boundary(const vec3 &pos) 
  {
    min = max = float4(pos.x, pos.y, pos.z, 0.0f);
  }
	
  const vec3 center() const {
		return float4(0.5f) * (max + min);
	}
	const vec3 hlen() const {
		return float4(0.5f) * (max - min);
	}

#ifdef __SSE_H__
	boundary(_v4sf _min, _v4sf _max)
  {
		min = _min;
		max = _max;
	}

	static const boundary merge(const boundary &a, const boundary &b){
		return boundary(
				__builtin_ia32_minps(a.min, b.min),
				__builtin_ia32_maxps(a.max, b.max));
	}
	void merge(const boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const boundary &a, const boundary &b)
  {
		return __builtin_ia32_movmskps(__builtin_ia32_orps(
         __builtin_ia32_cmpltps(a.max, b.min),
         __builtin_ia32_cmpltps(b.max, a.min)
         ));
	}
	friend bool overlapped(const boundary &a, const boundary &b){
		return !not_overlapped(a, b);
	}
#else
	static const boundary merge(const boundary &a, const boundary &b)
  {
		return boundary(mineach(a.min, b.min), maxeach(a.max, b.max));
	}
	void merge(const boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const boundary &a, const boundary &b){
		return (a.max.x() < b.min.x()) || (b.max.x() < a.min.x())
		    || (a.max.y() < b.min.y()) || (b.max.y() < a.min.y())
		    || (a.max.z() < b.min.z()) || (b.max.z() < a.min.z());
	}
	friend bool overlapped(const boundary &a, const boundary &b){
		return !not_overlapped(a, b);
	}
#endif
};
#if 0
template <>
struct boundary<float>{
	typedef vector3<float> vec;
	typedef float _v4sf __attribute__ ((vector_size(16)));

	vec min; 
	float p0;
	vec max;
	float p1;

	boundary() : min(HUGE), p0(0.f), max(-HUGE), p1(0.f) {}
	boundary(const vec &_min, const vec &_max) : min(_min), p0(0.f), max(_max), p1(0.f) {}
	boundary(const vec &pos, float h = 0.f) : 
		min(pos - vec(h)), p0(0.f), max(pos + vec(h)), p1(0.f) {}
	boundary(_v4sf _min, _v4sf _max){
		*(_v4sf *)&min = _min;
		*(_v4sf *)&max = _max;
	}

	static const boundary merge(const boundary &a, const boundary &b){
		return boundary(
				__builtin_ia32_minps(*(_v4sf *)&a.min, *(_v4sf *)&b.min),
				__builtin_ia32_maxps(*(_v4sf *)&a.max, *(_v4sf *)&b.max));
	}
	void merge(const boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const boundary &a, const boundary &b){
		return __builtin_ia32_movmskps(
				(_v4sf)(__builtin_ia32_cmpltps(
						*(_v4sf *)&a.max, *(_v4sf *)&b.min)))
		   ||  __builtin_ia32_movmskps(
				(_v4sf)(__builtin_ia32_cmpltps(
						*(_v4sf *)&b.max, *(_v4sf *)&a.min)));
	}
	friend bool overlapped(const boundary &a, const boundary &b){
		return !not_overlapped(a, b);
	}
};
#endif

#endif /* __BOUNDARY_H__ */
