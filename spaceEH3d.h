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

#ifndef SPACEEH3D_H
#define SPACEEH3D_H

#include <QtOpenGL>
#include <stdlib.h>

class Ccell3d;

class CSpaceEH3d
{
public:
  CSpaceEH3d(size_t sx, size_t sy, size_t sz);
  ~CSpaceEH3d();

  size_t sX, sY, sZ;  // size
  size_t sXY, sXYZ;  // size
  Ccell3d *c; // EH cells
  double *d;  // dither values
  GLuint YeeCell;

  void reset();
  void update_e();
  void update_h();
  void draw_yee();
};

#endif // SPACEEH3D_H
