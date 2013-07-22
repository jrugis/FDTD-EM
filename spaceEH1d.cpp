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

#include "defs.h"
#include "cell1d.h"
#include "spaceEH1d.h"

CSpaceEH1d::CSpaceEH1d(size_t s)
{
  size = s;
  c = new Ccell1d[size];  // EH cells
  eMax = new double[size];
  eMin = new double[size];
}

CSpaceEH1d::~CSpaceEH1d()
{
  delete[] c;
}

void CSpaceEH1d::reset()
{
  for(size_t i=0; i<size; i++) {
    c[i].reset();
    eMax[i] = 0.0;
    eMin[i] = 0.0;
  }
}

void CSpaceEH1d::update_e() // calculated using Ez Hy
{
#ifdef D22
  for (size_t i = 1; i < size; i++) {  // don't update lowest index e-field
    c[i].e = (c[i].cee * c[i].e) + (c[i].ceh * (c[i].h - c[i - 1].h));
    if(c[i].e > eMax[i]) eMax[i] = c[i].e;
    if(c[i].e < eMin[i]) eMin[i] = c[i].e;
  }
#endif
#ifdef D24
  for (size_t i = 2; i < size - 1; i++) {  // don't update lowest index e-field
    c[i].e = (c[i].cee * c[i].e) + (c[i].ceh * ((27 * (c[i].h - c[i-1].h)) + (c[i-2].h - c[i+1].h)) / 24);
  }
#endif
}

void CSpaceEH1d::update_h() // calculated using Ez Hy
{
#ifdef D22
  for (size_t i = 0; i < size - 1; i++) { // don't update highest index h-field
    c[i].h = (c[i].chh * c[i].h) + (c[i].che * (c[i + 1].e - c[i].e));
  }
#endif
#ifdef D24
    for (size_t i = 1; i < size - 2; i++) { // don't update highest index h-field
      c[i].h = (c[i].chh * c[i].h) + (c[i].che * ((27 * (c[i+1].e - c[i].e)) + (c[i-1].e - c[i+2].e)) / 24);
    }
#endif
}
