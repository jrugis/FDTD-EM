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

#include "spaceEH1d.h"
#include "cell1d.h"
#include "abc1o1d.h"

// abc constructed using e-field left, h-field right
CAbc1o1d::CAbc1o1d(CSpaceEH1d *s)
{
  size = s->size;
  c = s->c;
  double temp = sqrt(c[0].ceh * c[0].che);
  abcCoefL = (temp - 1.0) / (temp + 1.0);
  temp = sqrt(c[size - 1].ceh * c[size - 1].che);
  abcCoefR = (temp - 1.0) / (temp + 1.0);
}

void CAbc1o1d::reset()
{
  prevL = prevR = 0.0;
}

void CAbc1o1d::update_e() // left (e-field)
{
  c[0].e = prevL + abcCoefL * (c[1].e - c[0].e);
  prevL = c[1].e;
}
void CAbc1o1d::update_h() // right (h-field)
{
  c[size - 1].h = prevR + abcCoefR * (c[size - 2].h - c[size - 1].h);
  prevR = c[size - 2].h;
}
