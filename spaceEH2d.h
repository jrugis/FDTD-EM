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

#ifndef SPACEEH2D_H
#define SPACEEH2D_H

#include <stdlib.h>

class Ccell2d;

class CSpaceEH2d
{
public:
  CSpaceEH2d(size_t sx, size_t sy);
  ~CSpaceEH2d();

  size_t sX, sY;  // size
  size_t sXY;  // size
  Ccell2d *c; // EH cells
  double *d;  // dither values

  void reset();
  void update_e();
  void update_h();
};

#endif // SPACEEH2D_H
