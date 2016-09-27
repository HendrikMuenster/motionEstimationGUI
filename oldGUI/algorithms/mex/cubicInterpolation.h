#include "math.h"
#include "tools.h"
#include "linearInterpolation.h"

float cubicInterpolation(const float *image, float x, float y, const size_t *sizeMat);
float cubicInterpolate (float p[4], float x) ;
float bicubicInterpolate (float p[4][4], float x, float y);