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

#ifndef CELL2D_H
#define CELL2D_H

class Ccell2d
{
public:
  Ccell2d();
  void reset();

  double e, h1, h2;     // electric & magnetic fields
  double cee, ceh; // space parameters
  double ch1h, ch1e;
  double ch2h, ch2e;
};

#endif // CELL2D_H
