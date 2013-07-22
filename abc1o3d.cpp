/*
GL_10
An OpenGL+Qt4 FDTD electromagnetic simulation & visualization program.

Copyright (C) 2005-2012 John Rugis

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

rugis@msu.edu
*/

#include <math.h>

#include "defs.h"
#include "spaceEH3d.h"
#include "cell3d.h"
#include "abc1o3d.h"

// abc constructed using e-fields
CAbc1o3d::CAbc1o3d(CSpaceEH3d *s)
{
  sx = s->sX;
  sy = s->sY;
  sz = s->sZ;
  sxy = s->sXY;
  c = s->c;

  prevX0y = new double[sy * sz];
  prevX0z = new double[sy * sz];
  prevX1y = new double[sy * sz];
  prevX1z = new double[sy * sz];
  prevY0x = new double[sx * sz];
  prevY0z = new double[sx * sz];
  prevY1x = new double[sx * sz];
  prevY1z = new double[sx * sz];
  prevZ0x = new double[sx * sy];
  prevZ0y = new double[sx * sy];
  prevZ1x = new double[sx * sy];
  prevZ1y = new double[sx * sy];

  double temp = sqrt(c[0].cexh * c[0].chxe); // assumes uniform anisotropic
  abcCoef = (temp - 1.0) / (temp + 1.0);
}

void CAbc1o3d::reset()
{
  for(size_t i = 0; i < sy * sz; i++) {
    prevX0y[i] = prevX1y[i] = 0.0;
    prevX0z[i] = prevX1z[i] = 0.0;
  }
  for(size_t i = 0; i < sx * sz; i++) {
    prevY0x[i] = prevY1x[i] = 0.0;
    prevY0z[i] = prevY1z[i] = 0.0;
  }
  for(size_t i = 0; i < sx * sy; i++) {
    prevZ0x[i] = prevZ1x[i] = 0.0;
    prevZ0y[i] = prevZ1y[i] = 0.0;
  }
}

void CAbc1o3d::update_e()
{
  for (size_t k = 0; k < sz; k++) // ABC at "x0"
    for (size_t j = 0; j < sy; j++) {
      size_t m = j + k * sy;
      size_t n = j * sx + k * sxy;
      c[n].ey = prevX0y[m] + abcCoef * (c[n + 1].ey - c[n].ey);
      c[n].ez = prevX0z[m] + abcCoef * (c[n + 1].ez - c[n].ez);
      prevX0y[m] = c[n + 1].ey;
      prevX0z[m] = c[n + 1].ez;
    }
  for (size_t k = 0; k < sz; k++) // ABC at "x1"
    for (size_t j = 0; j < sy; j++) {
      size_t m = j + k * sy;
      size_t n = sx - 1 + j * sx + k * sxy;
      c[n].ey = prevX1y[m] + abcCoef * (c[n - 1].ey - c[n].ey);
      c[n].ez = prevX1z[m] + abcCoef * (c[n - 1].ez - c[n].ez);
      prevX1y[m] = c[n - 1].ey;
      prevX1z[m] = c[n - 1].ez;
    }
  for (size_t k = 0; k < sz; k++) // ABC at "y0"
    for (size_t i = 0; i < sx; i++) {
      size_t m = i + k * sx;
      size_t n = i + k * sxy;
      c[n].ex = prevY0x[m] + abcCoef * (c[n + sx].ex - c[n].ex);
      c[n].ez = prevY0z[m] + abcCoef * (c[n + sx].ez - c[n].ez);
      prevY0x[m] = c[n + sx].ex;
      prevY0z[m] = c[n + sx].ez;
    }
  for (size_t k = 0; k < sz; k++) // ABC at "y1"
    for (size_t i = 0; i < sx; i++) {
      size_t m = i + k * sx;
      size_t n = i + (sy - 1) * sx + k * sxy;
      c[n].ex = prevY1x[m] + abcCoef * (c[n - sx].ex - c[n].ex);
      c[n].ez = prevY1z[m] + abcCoef * (c[n - sx].ez - c[n].ez);
      prevY1x[m] = c[n - sx].ex;
      prevY1z[m] = c[n - sx].ez;
    }
  for (size_t j = 0; j < sy; j++) // ABC at "z0"
    for (size_t i = 0; i < sx; i++) {
      size_t m = i + j * sx;
      size_t n = i + j * sx;
      c[n].ex = prevZ0x[m] + abcCoef * (c[n + sxy].ex - c[n].ex);
      c[n].ey = prevZ0y[m] + abcCoef * (c[n + sxy].ey - c[n].ey);
      prevZ0x[m] = c[n + sxy].ex;
      prevZ0y[m] = c[n + sxy].ey;
    }
  for (size_t j = 0; j < sy; j++) // ABC at "z1"
    for (size_t i = 0; i < sx; i++) {
      size_t m = i + j * sx;
      size_t n = i + j * sx + (sz - 1) * sxy;
      c[n].ex = prevZ1x[m] + abcCoef * (c[n - sxy].ex - c[n].ex);
      c[n].ey = prevZ1y[m] + abcCoef * (c[n - sxy].ey - c[n].ey);
      prevZ1x[m] = c[n - sxy].ex;
      prevZ1y[m] = c[n - sxy].ey;
    }
}

void CAbc1o3d::update_h()
{
}
