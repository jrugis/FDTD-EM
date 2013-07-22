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

#include <stdlib.h>

#include "defs.h"
#include "utils.h"
#include "cell3d.h"
#include "spaceEH3d.h"

CSpaceEH3d::CSpaceEH3d(size_t sx, size_t sy, size_t sz)
{
  sX = sx;   // size
  sY = sy;
  sZ = sz;
  sXY = sx * sy;
  sXYZ = sx * sy * sz;
  c = new Ccell3d[sXYZ];  // EH cells
  d = new double[3 * sXYZ];   // dither values
  YeeCell = glGenLists(1);
  glNewList(YeeCell, GL_COMPILE);
    glLineWidth(2.0);
    glPointSize(6.0);
    glColor4d(MIDGREY, 0.5);
    glBegin(GL_LINE_LOOP);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.5, 0.0, 0.0);
      glVertex3d( 0.5, 0.5, 0.0);
      glVertex3d( 0.0, 0.5, 0.0);
    glEnd();
    glBegin(GL_LINE_LOOP);
      glVertex3d( 0.0, 0.0, 0.5);
      glVertex3d( 0.5, 0.0, 0.5);
      glVertex3d( 0.5, 0.5, 0.5);
      glVertex3d( 0.0, 0.5, 0.5);
    glEnd();
    glBegin(GL_LINES);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.0, 0.0, 0.5);
      glVertex3d( 0.5, 0.0, 0.0);
      glVertex3d( 0.5, 0.0, 0.5);
      glVertex3d( 0.5, 0.5, 0.0);
      glVertex3d( 0.5, 0.5, 0.5);
      glVertex3d( 0.0, 0.5, 0.0);
      glVertex3d( 0.0, 0.5, 0.5);
    glEnd();
    glColor4d(MIDYELLOW, 1.0);
    glBegin(GL_LINES);
      glVertex3d( 0.45, 0.0, 0.0);
      glVertex3d( 0.55, 0.0, 0.0);
      glVertex3d( 0.0, 0.45, 0.0);
      glVertex3d( 0.0, 0.55, 0.0);
      glVertex3d( 0.0, 0.0, 0.45);
      glVertex3d( 0.0, 0.0, 0.55);
    glEnd();
    glBegin(GL_POINTS);
      glVertex3d( 0.45, 0.0, 0.0);
      glVertex3d( 0.0, 0.45, 0.0);
      glVertex3d( 0.0, 0.0, 0.45);
    glEnd();
    glColor4d(MIDCYAN, 1.0);
    glBegin(GL_LINES);
      glVertex3d(-0.05, 0.5, 0.5);
      glVertex3d(0.05, 0.5, 0.5);
      glVertex3d(0.5, -0.05, 0.5);
      glVertex3d(0.5, 0.05, 0.5);
      glVertex3d(0.5, 0.5, -0.05);
      glVertex3d(0.5, 0.5, 0.05);
    glEnd();
    glBegin(GL_POINTS);
      glVertex3d(-0.05, 0.5, 0.5);
      glVertex3d(0.5, -0.05, 0.5);
      glVertex3d(0.5, 0.5, -0.05);
    glEnd();
  glEndList();
}

CSpaceEH3d::~CSpaceEH3d()
{
    delete[] c;
}

void CSpaceEH3d::reset()
{
  for(size_t i = 0; i < sXYZ; i++) c[i].reset();
  for(size_t i = 0; i < 3 * sXYZ; i++) d[i] = randpm();
}

void CSpaceEH3d::draw_yee()
{
  glCallList(YeeCell);
}


void CSpaceEH3d::update_e()
{
  for (size_t k = 1; k < sZ - 1; k++) { // don't update boundary e-fields
    for (size_t j = 1; j < sY - 1; j++) {
      for (size_t i = 1; i < sX - 1; i++) {
        size_t n = i + j * sX + k * sXY;
        c[n].ex =
            c[n].cexe * c[n].ex
          + c[n].cexh * ((c[n].hz - c[n - sX].hz) - (c[n].hy - c[n - sXY].hy));
        c[n].ey =
            c[n].ceye * c[n].ey
          + c[n].ceyh * ((c[n].hx - c[n - sXY].hx) - (c[n].hz - c[n - 1].hz));
        c[n].ez =
            c[n].ceze * c[n].ez
          + c[n].cezh * ((c[n].hy - c[n - 1].hy) - (c[n].hx - c[n - sX].hx));
      }
    }
  }
}

void CSpaceEH3d::update_h()
{
  for (size_t k = 0; k < sZ - 1; k++) {  // don't update highest h-fields
    for (size_t j = 0; j < sY - 1; j++) {
      for (size_t i = 0; i < sX - 1; i++) {
        size_t n = i + j * sX + k * sXY;
        c[n].hx =
            c[n].chxh * c[n].hx
          + c[n].chxe * ((c[n + sXY].ey - c[n].ey) - (c[n + sX].ez - c[n].ez));
        c[n].hy =
            c[n].chyh * c[n].hy
          + c[n].chye * ((c[n + 1].ez - c[n].ez) - (c[n + sXY].ex - c[n].ex));
        c[n].hz =
            c[n].chzh * c[n].hz
          + c[n].chze * ((c[n + sX].ex - c[n].ex) - (c[n + 1].ey - c[n].ey));
      }
    }
  }
}

