#ifndef __PRIMITIVES_H__
#define __PRIMITIVES_H__

template<class T>
static inline T __min(const T &a, const T &b) {return a < b ? a : b;}
template<class T>
static inline T __max(const T &a, const T &b) {return a > b ? a : b;}

#define SQR(x) ((x)*(x))

#include "vector3.h"

typedef float real;
typedef vector3<real> vec3;

#endif /* __PRIMITIVES_H__ */
