#ifndef _OPENMHD_PARAM_
#define _OPENMHD_PARAM_
#include <math.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define INDEX(j,i) ((j)*ix+(i))

namespace param {
  const int version = 20191227;
  const double gamma = 5.0 / 3.0;
  const double pi = 4.0 * atan(1.0);

  const int var1 = 9, var2 = 4;
  const int mx = 0, vx = 0;
  const int my = 1, vy = 1;
  const int mz = 2, vz = 2;
  const int en = 3, pr = 3;
  const int ro = 4;
  const int bx = 5;
  const int by = 6;
  const int bz = 7;
  const int ps = 8;
}

# endif
