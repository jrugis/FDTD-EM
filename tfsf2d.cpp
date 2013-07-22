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

#include "defs.h"
#include "cell1d.h"
#include "cell2d.h"
#include "spaceEH1d.h"
#include "spaceEH2d.h"
#include "abc2o1d.h"
#include "tfsf2d.h"

CTfsf2d::CTfsf2d(CSpaceEH2d *s, size_t sB, size_t sD)
{
  #define MAX_LOSS 0.35
  c = s->c;
  sx = s->sX;
  sy = s->sY;
  sb = sB;   // tfsf boundary (per edge)
  sd = sD;   // tfsf aux start & decay region
  sa = 2 * sd + sx;
  a = new CSpaceEH1d(sa);
  inp = &(a->c[sd].e);       // source input
  inpm1 = &(a->c[sd - 1].h); // tfsf

  // setup aux material
  for (size_t i = 0; i < sx; i++) { // copy material strip from 2d model
    size_t m = sd + i;
    a->c[m].cee = c[i].cee;
    a->c[m].ceh = c[i].ceh;
    a->c[m].chh = c[i].ch2h;  // use Hy values
    a->c[m].che = c[i].ch2e;
  }
  for (size_t i = 0; i < sd; i++) { // LHS duplicate material
    a->c[i].cee = a->c[sd].cee;
    a->c[i].ceh = a->c[sd].ceh;
    a->c[i].chh = a->c[sd].chh;
    a->c[i].che = a->c[sd].che;
  }
  for (size_t i = 0; i < sd; i++) { // RHS duplicate material & smooth loss
    size_t m = sd + sx + i;
    size_t n = sd + sx - 1;
    double lossFactor = MAX_LOSS * pow((i + 0.5) / sd, 2);  // fractional depth squared
    a->c[m].cee = a->c[n].cee * (1.0 - lossFactor) / (1.0 + lossFactor);
    a->c[m].ceh = a->c[n].ceh / (1.0 + lossFactor);
    lossFactor = MAX_LOSS * pow((i + 1.0) / sd, 2); // h field is offset (deeper) by 0.5
    a->c[m].chh = a->c[n].chh * (1.0 - lossFactor) / (1.0 + lossFactor);
    a->c[m].che = a->c[n].che / (1.0 + lossFactor);
  }
  // set abc's after material initialization!!!
  abc2o1d = new CAbc2o1d(a);
}

void CTfsf2d::reset()
{
  a->reset();
  abc2o1d->reset();
}

void CTfsf2d::updateA()
{
  // correct Hy along left edge
  size_t i = sb - 1;
  for (size_t j = sb; j < sy - sb; j++) {
    size_t n = i + j * sx;
    c[n].h2 -= c[n].ch2e * a->c[sd + i + 1].e; // h2 is y direction
  }
  // correct Hy along right edge
  i = sx - sb - 1;
  for (size_t j = sb; j < sy - sb; j++) {
    size_t n = i + j * sx;
    c[n].h2 += c[n].ch2e * a->c[sd + i].e; // h2 is y direction
  }
  // correct Hx along the bottom
  size_t j = sb - 1;
  for (size_t i = sb; i < sx - sb; i++) {
    size_t n = i + j * sx;
    c[n].h1 += c[n].ch1e * a->c[sd + i].e; // h1 is x direction
  }
  // correct Hx along the top
  j = sy - sb - 1;
  for (size_t i = sb; i < sx - sb; i++) {
    size_t n = i + j * sx;
    c[n].h1 -= c[n].ch1e * a->c[sd + i].e; // h1 is x direction
  }

  a->update_h(); // update magnetic field
}

void CTfsf2d::updateB()
{
  abc2o1d->update_h();
  a->update_e(); // update electric field
  abc2o1d->update_e();

  // correct Ez field along left edge
  size_t i = sb;
  for (size_t j = sb; j < sy - sb; j++) {
    size_t n = i + j * sx;
    c[n].e -= c[n].ceh * a->c[sd + i - 1].h;
  }
  // correct Ez field along right edge
  i = sx - sb - 1;
  for (size_t j = sb; j < sy - sb; j++) {
    size_t n = i + j * sx;
    c[n].e += c[n].ceh * a->c[sd + i].h;
  }
}
