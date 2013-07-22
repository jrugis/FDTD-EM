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

#include "utils.h"
#include "cell2d.h"
#include "spaceEH2d.h"

CSpaceEH2d::CSpaceEH2d(size_t sx, size_t sy)
{
  sX = sx;   // size
  sY = sy;
  sXY = sx * sy;
  c = new Ccell2d[sXY];  // EH cells
  d = new double[3 * sXY];   // dither values

}

CSpaceEH2d::~CSpaceEH2d()
{
  delete[] c;
}

void CSpaceEH2d::reset()
{
  for(size_t i = 0; i < sXY; i++) c[i].reset();
  for(size_t i = 0; i < 3 * sXY; i++) d[i] = randpm();
}

void CSpaceEH2d::update_e() // calculated using Ez Hx Hy
{
  for (size_t j = 1; j < sY; j++) { // don't update lowest e-fields
    for (size_t i = 1; i < sX; i++) {
      size_t n = i + j * sX;
      c[n].e =
          c[n].cee * c[n].e
        + c[n].ceh * ((c[n].h2 - c[n - 1].h2) - (c[n].h1 - c[n - sX].h1));
    }
  }
}

void CSpaceEH2d::update_h() // calculated using Ez Hx Hy
{
  for (size_t j = 0; j < sY - 1; j++) {  // don't update highest h-fields
    for (size_t i = 0; i < sX - 1; i++) {
      size_t n = i + j * sX;
      c[n].h1 =
          c[n].ch1h * c[n].h1
        - c[n].ch1e * (c[n + sX].e - c[n].e);
      c[n].h2 =
          c[n].ch2h * c[n].h2
        + c[n].ch2e * (c[n + 1].e - c[n].e);
    }
  }
}
