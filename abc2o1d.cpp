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
#include "abc2o1d.h"

CAbc2o1d::CAbc2o1d(CSpaceEH1d *s)
{
  size = s->size;
  c = s->c;

  // left coefficients
  double temp1 = sqrt(s->c[0].ceh * s->c[0].che);
  double temp2 = 1.0 / temp1 + 2.0 + temp1;
  abcCoefL[0] = -(1.0 / temp1 - 2.0 + temp1) / temp2;
  abcCoefL[1] = -2.0 * (temp1 - 1.0 / temp1) / temp2;
  abcCoefL[2] = 4.0 * (temp1 + 1.0 / temp1) / temp2;

  // right coefficients
  temp1 = sqrt(s->c[(s->size) - 1].ceh * s->c[(s->size) - 1].che);
  temp2 = 1.0 / temp1 + 2.0 + temp1;
  abcCoefR[0] = -(1.0 / temp1 - 2.0 + temp1) / temp2;
  abcCoefR[1] = -2.0 * (temp1 - 1.0 / temp1) / temp2;
  abcCoefR[2] = 4.0 * (temp1 + 1.0 / temp1) / temp2;
}

void CAbc2o1d::reset()
{
    for (int j = 0; j < 2; j++) // time: back
      for (int i = 0; i < 3; i++) // space: from edge
          prevL[i][j] = prevR[i][j] = 0.0;
}

void CAbc2o1d::update_e()  // left (e-field)
{
  c[0].e =
      abcCoefL[0] * (c[2].e + prevL[0][1])
    + abcCoefL[1] * (prevL[0][0] + prevL[2][0] - c[1].e - prevL[1][1])
    + abcCoefL[2] * prevL[1][0] - prevL[2][1];

  for (int i = 0; i < 3; i++) {
    prevL[i][1] = prevL[i][0];
    prevL[i][0] = c[i].e;
  }
}
void CAbc2o1d::update_h()  // right (h-field)
{
  c[size - 1].h =
      abcCoefR[0] * (c[size - 3].h + prevR[0][1])
    + abcCoefR[1] * (prevR[0][0] + prevR[2][0] - c[size - 2].h - prevR[1][1])
    + abcCoefR[2] * prevR[1][0] - prevR[2][1];

  for (int i = 0; i < 3; i++) {
    prevR[i][1] = prevR[i][0];
    prevR[i][0] = c[size - 1 - i].h;
  }
}
