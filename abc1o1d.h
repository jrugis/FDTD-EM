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

#ifndef ABC1O1D_H
#define ABC1O1D_H

class Ccell1d;
class CSpaceEH1d;

class CAbc1o1d
{
public:
  CAbc1o1d(CSpaceEH1d *s);

  void reset();
  void update_e();
  void update_h();

private:
  size_t size;
  Ccell1d *c;
  double prevL, prevR;
  double abcCoefL, abcCoefR;
};

#endif // ABC1O1D_H
