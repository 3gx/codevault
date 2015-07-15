#include <sys/time.h>

#ifndef __TIMER__H__
#define __TIMER__H__

static double get_wtime(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#endif // __TIMER__H__
