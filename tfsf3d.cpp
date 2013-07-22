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
#include "cell3d.h"
#include "spaceEH1d.h"
#include "spaceEH3d.h"
#include "abc2o1d.h"
#include "tfsf3d.h"

CTfsf3d::CTfsf3d(CSpaceEH3d *s, size_t sB, size_t sD)
{
  #define MAX_LOSS 0.35
  c = s->c;
  sx = s->sX;
  sy = s->sY;
  sz = s->sZ;
  sxy = sx * sy;
  sb = sB;   // tfsf boundary (per edge)
  sd = sD;   // tfsf aux start & decay region
//  sa = 2 * sd + sx;
  sa = 2 * sd + sz;
  a = new CSpaceEH1d(sa);
  inp = &(a->c[sd].e);       // source input
  inpm1 = &(a->c[sd - 1].h); // tfsf

  // setup aux material
//  for (size_t i = 0; i < sx; i++) { // copy material strip from 3d model
  for (size_t i = 0; i < sz; i++) { // copy material strip from 3d model
    size_t m = sd + i;
    size_t n = i * sxy;
//    a->c[m].cee = c[i].ceze;
//    a->c[m].ceh = c[i].cezh;
//    a->c[m].chh = c[i].chyh;
//    a->c[m].che = c[i].chye;
    a->c[m].cee = c[n].ceze;
    a->c[m].ceh = c[n].cezh;
    a->c[m].chh = c[n].chyh;
    a->c[m].che = c[n].chye;

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

void CTfsf3d::reset()
{
  a->reset();
  abc2o1d->reset();
}

void CTfsf3d::updateA()
{
  // correct Hy at low x
  size_t i = sb;
  for (size_t k = sb; k < sz - sb; k++) {
    for (size_t j = sb; j <= sy - sb; j++) {
      size_t n = i + j * sx + k * sxy;
      c[n - 1].hy -= c[n].chye * a->c[sd + i].e;
//      c[n - 1].hy -= c[n].chye * a->c[sd + k].e;
    }
  }
/*// correct Hy at firstX-1/2 by subtracting Ez_inc
mm = firstX;
for (nn = firstY; nn <= lastY; nn++)
for (pp = firstZ; pp < lastZ; pp++
Hy(mm - 1, nn, pp) -= Chye(mm, nn, pp) * Ez1G(g1, mm);*/

  // correct Hy at high x
  i = sx - sb;
  for (size_t k = sb; k < sz - sb; k++) {
    for (size_t j = sb; j <= sy - sb; j++) {
      size_t n = i + j * sx + k * sxy;
      c[n].hy += c[n].chye * a->c[sd + i].e;
//      c[n].hy += c[n].chye * a->c[sd + k].e;
    }
  }
/*// correct Hy at lastX + 1/2 by adding Ez_inc
mm = lastX;
for (nn = firstY; nn <= lastY; nn++)
for (pp = firstZ; pp < lastZ; pp++)
Hy(mm, nn, pp) += Chye(mm, nn, pp) * Ez1G(g1, mm);*/

  // correct Hx at low y
  size_t j = sb - 1;
  for (size_t k = sb; k < sz - sb; k++) {
    for (size_t i = sb; i <= sx - sb; i++) {
      size_t n = i + j * sx + k * sxy;
      c[n].hx += c[n].chxe * a->c[sd + i].e;
//      c[n].hx += c[n].chxe * a->c[sd + k].e;
    }
  }
// correct Hx at firstY-1/2 by adding Ez_inc
/*nn = firstY;
for (mm = firstX; mm <= lastX; mm++)
for (pp = firstZ; pp < lastZ; pp++)
Hx(mm, nn - 1, pp) += Chxe(mm, nn - 1, pp) * Ez1G(g1, mm);*/

  // correct Hx at high y
  j = sy - sb;
  for (size_t k = sb; k < sz - sb; k++) {
    for (size_t i = sb; i <= sx - sb; i++) {
      size_t n = i + j * sx + k * sxy;
      c[n].hx -= c[n].chxe * a->c[sd + i].e;
//      c[n].hx -= c[n].chxe * a->c[sd + k].e;
    }
  }
/*// correct Hx at lastY+1/2 by subtracting Ez_inc
nn = lastY;
for (mm = firstX; mm <= lastX; mm++)
for (pp = firstZ; pp < lastZ; pp++)
Hx(mm, nn, pp) -= Chxe(mm, nn, pp) * Ez1G(g1,mm);*/

  a->update_h(); // update magnetic field
}

void CTfsf3d::updateB()
{
  abc2o1d->update_h();
  a->update_e(); // update electric field
  abc2o1d->update_e();

  // correct Ez field at low x
  size_t i = sb;
  for (size_t k = sb; k < sz - sb; k++) {
    for (size_t j = sb; j <= sy - sb; j++) {
      size_t n = i + j * sx + k * sxy;
      c[n].ez -= c[n].cezh * a->c[sd + i - 1].h;
//      c[n].ez -= c[n].cezh * a->c[sd + k - 1].h;
    }
  }
// correct Ez at firstX face by subtracting Hy_inc
/*mm = firstX;
for (nn = firstY; nn <= lastY; nn++)
for (pp = firstZ; pp < lastZ; pp++)
Ez(mm, nn, pp) -= Cezh(mm, nn, pp) * Hy1G(g1, mm - 1);*/

  // correct Ez field at high x
  i = sx - sb;
  for (size_t k = sb; k < sz - sb; k++) {
    for (size_t j = sb; j <= sy - sb; j++) {
      size_t n = i + j * sx + k * sxy;
      c[n].ez += c[n].cezh * a->c[sd + i].h;
//      c[n].ez += c[n].cezh * a->c[sd + k].h;
    }
  }
// correct Ez at lastX face by adding Hy_inc
/*mm = lastX;
for (nn = firstY; nn <= lastY; nn++)
for (pp = firstZ; pp < lastZ; pp++)
Ez(mm, nn, pp) += Cezh(mm, nn, pp) * Hy1G(g1, mm);*/

  // correct Ex field at low z
  size_t k = sb;
  for (size_t j = sb; j <= sy - sb; j++) {
    for (size_t i = sb; i < sx - sb; i++) {
      size_t n = i + j * sx + k * sxy;
      c[n].ex += c[n].cexh * a->c[sd + i].h;
//      c[n].ex += c[n].cexh * a->c[sd + k].h;
    }
  }
// correct Ex at firstZ face by adding Hy_inc
/*pp = firstZ;
for (mm = firstX; mm < lastX; mm++)
for (nn = firstY; nn <= lastY; nn++)
Ex(mm, nn, pp) += Cexh(mm, nn, pp) * Hy1G(g1, mm);*/

  // correct Ex field at high z
  k = sz - sb;
  for (size_t j = sb; j <= sy - sb; j++) {
    for (size_t i = sb; i < sx - sb; i++) {
      size_t n = i + j * sx + k * sxy;
      c[n].ex -= c[n].cexh * a->c[sd + i].h;
//      c[n].ex -= c[n].cexh * a->c[sd + k].h;
    }
  }
/*// correct Ex at lastZ face by subtracting Hy_inc
pp = lastZ;
for (mm = firstX; mm < lastX; mm++)
for (nn = firstY; nn <= lastY; nn++)
Ex(mm, nn, pp) -= Cexh(mm, nn, pp) * Hy1G(g1, mm);*/
}
