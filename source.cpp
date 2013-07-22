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
#include "source.h"
#include "spaceEH1d.h"

void sourceGaussian(bool a, double s, double *eh, double ts, double dts, double nwtss)
{
  double temp = (a ? *eh : 0.0); // additive?
  *eh = temp + s * exp( pow(ts - dts, 2) / nwtss);
}

void sourceSine(bool a, double s, double *eh, size_t ts, double omega)
{
  double temp = (a ? *eh : 0.0); // additive?
  *eh = temp + s * sin(omega * ts);
}

void sourceRicker(bool a, double s, double *eh, size_t ts, double wts)
{
  double temp = (a ? *eh : 0.0); // additive?
  double arg = pow(2* PI * (ts / wts - 1.0), 2);
  *eh = temp + s * (1.0 - 2.0 * arg) * exp(-arg);
}

void sourceImpulse(bool a, double s, double *eh, size_t ts, size_t on, size_t off)
{
  double temp = (a ? *eh : 0.0); // additive?
  *eh = temp + ((ts >= on && ts < off) ? s : 0.0);
}
