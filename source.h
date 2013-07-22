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

#ifndef SOURCE_H
#define SOURCE_H

#include <stdlib.h>

/*
  signal source arguments:
   1 additive
   2 scaling factor
   3 destination pointer
   4 time step
   5 other... (source specific)
*/

void sourceGaussian(bool additive, double scaling, double *eh, double time_step, double dts, double nwtss);
void sourceSine(bool additive, double scaling, double *eh, size_t time_step, double omega);
void sourceRicker(bool additive, double scaling, double *eh, size_t time_step, double wts);
void sourceImpulse(bool additive, double scaling, double *eh, size_t time_step, size_t on, size_t off);
#endif // SOURCE_H
