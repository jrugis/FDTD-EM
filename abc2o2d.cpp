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

#include "cell2d.h"
#include "spaceEH2d.h"
#include "abc2o2d.h"

CAbc2o2d::CAbc2o2d(CSpaceEH2d *s)
{
  sX = s->sX;
  sY = s->sY;
  c = s->c;
  for (int j = 0; j < 2; j++) // time: back
    for (int i = 0; i < 3; i++) { // position: from edge
      prevL[i][j] = new double[sY];
      prevR[i][j] = new double[sY];
      prevT[i][j] = new double[sX];
      prevB[i][j] = new double[sX];
    }

  // values from one corner used, assumes homogeneous boundary
  double temp1 = sqrt(c[0].ceh * c[0].ch1e);
  double temp2 = 1.0 / temp1 + 2.0 + temp1;
  coef[0] = -(1.0 / temp1 - 2.0 + temp1) / temp2;
  coef[1] = -2.0 * (temp1 - 1.0 / temp1) / temp2;
  coef[2] = 4.0 * (temp1 + 1.0 / temp1) / temp2;
}

void CAbc2o2d::reset()
{
  for (unsigned int j = 0; j < 2; j++) // time: back
    for (unsigned int i = 0; i < 3; i++) { // position: from edge
      for (unsigned int k = 0; k < sY; k++)
        prevL[i][j][k] = prevR[i][j][k] = 0.0;
      for (unsigned int k = 0; k < sX; k++)
        prevT[i][j][k] = prevB[i][j][k] = 0.0;
    }
}

void CAbc2o2d::update()
{
  for (unsigned int k = 0; k < sY; k++) {  // left
    c[k * sX].e =
        coef[0] * (c[2 + k * sX].e + prevL[0][1][k])
      + coef[1] * (prevL[0][0][k] + prevL[2][0][k] - c[1 + k * sX].e - prevL[1][1][k])
      + coef[2] * prevL[1][0][k] - prevL[2][1][k];
    for (unsigned int i = 0; i < 3; i++) { // store prev fields
      prevL[i][1][k] = prevL[i][0][k];
      prevL[i][0][k] = c[i + k * sX].e;
    }
  }
  for (unsigned int k = 0; k < sY; k++) {  // right
    c[sX - 1 + k * sX].e =
        coef[0] * (c[sX - 3 + k * sX].e + prevR[0][1][k])
      + coef[1] * (prevR[0][0][k] + prevR[2][0][k] - c[sX - 2 + k * sX].e - prevR[1][1][k])
      + coef[2] * prevR[1][0][k] - prevR[2][1][k];
    for (unsigned int i = 0; i < 3; i++) { // store prev fields
      prevR[i][1][k] = prevR[i][0][k];
      prevR[i][0][k] = c[sX - 1 - i + k * sX].e;
    }
  }
  for (unsigned int k = 0; k < sX; k++) {  // bottom
    c[k].e =
        coef[0] * (c[k + 2 * sX].e + prevB[0][1][k])
      + coef[1] * (prevB[0][0][k] + prevB[2][0][k] - c[k + sX].e - prevB[1][1][k])
      + coef[2] * prevB[1][0][k] - prevB[2][1][k];
    for (unsigned int i = 0; i < 3; i++) { // store prev fields
      prevB[i][1][k] = prevB[i][0][k];
      prevB[i][0][k] = c[k + i * sX].e;
    }
  }
  for (unsigned int k = 0; k < sX; k++) {  // top
    c[k + (sY - 1) * sX].e =
        coef[0] * (c[k + (sY - 3) * sX].e + prevT[0][1][k])
      + coef[1] * (prevT[0][0][k] + prevT[2][0][k] - c[k + (sY - 2) * sX].e - prevT[1][1][k])
      + coef[2] * prevT[1][0][k] - prevT[2][1][k];
    for (unsigned int i = 0; i < 3; i++) { // store prev fields
      prevT[i][1][k] = prevT[i][0][k];
      prevT[i][0][k] = c[k + (sY - 1 - i) * sX].e;
    }
  }
}

