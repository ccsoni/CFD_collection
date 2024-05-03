#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <float.h>

typedef uint64_t UINT;
typedef double REAL;

#define PI (3.14159265359)

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

struct run_param {
  char model_name[32];
  char scheme_name[32];

  UINT step;
  UINT nmesh;
  REAL xmin, xmax, delta_x;
  REAL tnow, tend, dtime;

  REAL CFL;
};
