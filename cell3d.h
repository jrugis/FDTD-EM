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

#ifndef CELL3D_H
#define CELL3D_H

class Ccell3d
{
public:
  Ccell3d();
  void reset();

  double ex, ey, ez, hx, hy, hz;     // electric & magnetic fields
  double cexe, cexh; // space parameters
  double ceye, ceyh;
  double ceze, cezh;
  double chxh, chxe;
  double chyh, chye;
  double chzh, chze;
};

#endif // CELL3D_H
