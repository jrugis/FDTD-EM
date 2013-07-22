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

#include <fstream>
#include <limits.h>

#include "defs.h"
#include "utils.h"
#include "glwidget.h"
#include "source.h"
#include "cell2d.h"
#include "spaceEH2d.h"
#include "abc2o2d.h"
#include "tfsf2d.h"

#include "model2d.h"

#define SIZEX 101
#define SIZEY 101
#define SIZEZ 101

#define FREE_SPACE
//#define HALF_SPACE_LOSSLESS_DIELECTRIC

#define PEC_DISK_00
//#define PEC_DISK_01
//#define PEC_BOX
//#define PEC_LINE
//#define PEC_SLIT

#define RICKER_PLANE
//#define RICKER_POINT
//#define GAUSSIAN_POINT
//#define SNAPSHOT

#define ABC_SECOND_ORDER

CModel2D::CModel2D(GLWidget *parent) : CModel(parent)
{
  space2d = NULL;
  abc2o2d = NULL;
  tfsf2d = NULL;  // tfsf's

  set_material(); // space & material
  reset();        // space & material
}

// ***********************************************************************
// model materials
// ***********************************************************************
void CModel2D::set_material()
{
  space2d = new CSpaceEH2d(SIZEX, SIZEY);

#if defined FREE_SPACE
  for(size_t j = 0; j < SIZEY; j++) {
    for(size_t i = 0; i < SIZEX; i++) {
      size_t m = i + j * SIZEX;
      space2d->c[m].cee = 1.0; // free-space
      space2d->c[m].ceh = DTDS2D * IMP0;
      space2d->c[m].ch1h = 1.0;
      space2d->c[m].ch1e = DTDS2D / IMP0;
      space2d->c[m].ch2h = 1.0;
      space2d->c[m].ch2e = DTDS2D / IMP0;
    }
  }
#endif

#if defined HALF_SPACE_LOSSLESS_DIELECTRIC
  #define RHS_EPSR 8.0
  for(size_t j = 0; j < SIZEY; j++) {
    for(size_t i = 0; i < SIZEX / 2; i++) {
      size_t m = i + j * SIZEX;
      space2d->c[m].cee = 1.0; // free-space
      space2d->c[m].ceh = DTDS2D * IMP0;
      space2d->c[m].ch1h = 1.0;
      space2d->c[m].ch1e = DTDS2D / IMP0;
      space2d->c[m].ch2h = 1.0;
      space2d->c[m].ch2e = DTDS2D / IMP0;
    }
  }
  for(size_t j = 0; j < SIZEY; j++) {
    for(size_t i = SIZEX / 2; i < SIZEX; i++) {
      size_t m = i + j * SIZEX;
      space2d->c[m].cee = 1.0;
      space2d->c[m].ceh = DTDS2D * IMP0 / RHS_EPSR;
      space2d->c[m].ch1h = 1.0;
      space2d->c[m].ch1e = DTDS2D / IMP0;
      space2d->c[m].ch2h = 1.0;
      space2d->c[m].ch2e = DTDS2D / IMP0;
    }
  }
#endif

#ifdef PEC_DISK_00
  #define CR (0.15 * SIZEY)
  #define CX (0.5 * SIZEX)
  #define CY (0.5 * SIZEY)
  for(int j = CY - CR; j < CY + CR; j++) {
    for(int i = CX - CR; i < CX + CR; i++) {
      int cr2 = CR * CR;
      if( pow((i - CX), 2) + pow((j - CY), 2) <= cr2) {
        int n = i + j * SIZEX;
        space2d->c[n].cee = 0.0;
        space2d->c[n].ceh = 0.0;
      }
    }
  }
#endif

#ifdef PEC_DISK_01
  #define CR (0.05 * SIZEX)
  #define CX (0.75 * SIZEX)
  #define CY (0.25 * SIZEY)
  for(int j = CY - CR; j < CY + CR; j++) {
    for(int i = CX - CR; i < CX + CR; i++) {
      int cr2 = CR * CR;
      if( pow((i - CX), 2) + pow((j - CY), 2) <= cr2) {
        int n = i + j * SIZEX;
        space2d->c[n].cee = 0.0;
        space2d->c[n].ceh = 0.0;
      }
    }
  }
#endif

#ifdef PEC_BOX
  for(int j = 0; j < 200; j++) {
    space2d->c[250 + (j +150) * SIZEX].cee = 0.0;
    space2d->c[250 + (j +150) * SIZEX].ceh = 0.0;
    space2d->c[350 + (j +150) * SIZEX].cee = 0.0;
    space2d->c[350 + (j +150) * SIZEX].ceh = 0.0;
    //space2d->c[(j + 250) + 150 * SIZEX].cee = 0.0;
    //space2d->c[(j + 250) + 150 * SIZEX].ceh = 0.0;
    //space2d->c[(j + 250) + 250 * SIZEX].cee = 0.0;
    //space2d->c[(j + 250) + 250 * SIZEX].ceh = 0.0;
  }
#endif

#ifdef PEC_LINE
  #define CS (SIZEY / 2)
  #define CX (SIZEX / 3)
  #define CY1 (SIZEY / 2 - CS / 2)
  #define CY2 (CY1 + CS)
  for(int j = CY1; j < CY2; j++) {
    space2d->c[CX + j * SIZEX].cee = 0.0;
    space2d->c[CX + j * SIZEX].ceh = 0.0;
  }
#endif

#ifdef PEC_SLIT
  #define CS 6
  #define CX (SIZEX / 3)
  #define CY1 (SIZEY / 2 - CS / 2)
  #define CY2 (CY1 + CS)
  for(int j = 0; j < CY1; j++) {
    space2d->c[CX + j * SIZEX].cee = 0.0;
    space2d->c[CX + j * SIZEX].ceh = 0.0;
  }
  for(int j = CY2; j < SIZEY; j++) {
    space2d->c[CX + j * SIZEX].cee = 0.0;
    space2d->c[CX + j * SIZEX].ceh = 0.0;
  }
#endif

  // set tfsf and abc's after material initialization!!!

#ifdef RICKER_PLANE
  #define SD 10  // tfsf aux decay size
  #define SB 3  // tfsf boundary size
  tfsf2d = new CTfsf2d(space2d, SB, SD);
#endif

#ifdef ABC_SECOND_ORDER
  abc2o2d = new CAbc2o2d(space2d);
#endif

}

// ***********************************************************************
// model field step
// ***********************************************************************
void CModel2D::step()
{
  space2d->update_h(); // ***** update magnetic field *****

#ifdef RICKER_PLANE
  #define WTS 50
  tfsf2d->updateA();
  sourceRicker(true, 0.3 / -IMP0, tfsf2d->inpm1, time_step, WTS);
  sourceRicker(true, 0.3, tfsf2d->inp, time_step + 1, WTS);
  tfsf2d->updateB();
#endif

  space2d->update_e();  // ***** update electric field *****

#ifdef RICKER_POINT
  #define WTS 200
  #define SRC_X (SIZEX / 4)
  #define SRC_Y (SIZEY / 2)
  #define SRC (SRC_X + SRC_Y * SIZEX)
  sourceRicker(false, 0.8, &(space2d->c[SRC].e), time_step, WTS);
#endif

#ifdef GAUSSIAN_POINT
  #define WTS 40
  #define DTS (WTS * 4) // delay time steps
  #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  #define AMPE  0.5
  #define AMPH  (AMPE / IMP0)
  #define SRC_X (SIZEX / 2)
  #define SRC_Y (SIZEY / 2)
  #define SRC0 (SRC_X +      SRC_Y * SIZEX)
  #define SRC1 (SRC_X + 1 +  SRC_Y * SIZEX)
  #define SRC2 (SRC_X + 2 +  SRC_Y * SIZEX)
  #define SRC3 (SRC_X + 3 +  SRC_Y * SIZEX)
  #define SRC4 (SRC_X +     (SRC_Y + 1) * SIZEX)
  #define SRC5 (SRC_X + 1 + (SRC_Y + 1) * SIZEX)
  #define SRC6 (SRC_X + 2 + (SRC_Y + 1) * SIZEX)
  #define SRC7 (SRC_X + 3 + (SRC_Y + 1) * SIZEX)
  #define SRC8 (SRC_X +     (SRC_Y + 2) * SIZEX)
  #define SRC9 (SRC_X + 1 + (SRC_Y + 2) * SIZEX)
  #define SRCA (SRC_X + 2 + (SRC_Y + 2) * SIZEX)
  #define SRCB (SRC_X + 3 + (SRC_Y + 2) * SIZEX)

  sourceGaussian(true, AMPE, &(space2d->c[SRC0].e), time_step, DTS, NWTSS);
  sourceGaussian(true, AMPE, &(space2d->c[SRC1].e), time_step, DTS, NWTSS);
  sourceGaussian(true, AMPE, &(space2d->c[SRC4].e), time_step, DTS, NWTSS);
  sourceGaussian(true, AMPE, &(space2d->c[SRC5].e), time_step, DTS, NWTSS);

//  sourceGaussian(true, -AMPH, &(space2d->c[SRC1].h2), time_step, DTS, NWTSS);
//  sourceGaussian(true, -AMPH, &(space2d->c[SRC4].h1), time_step, DTS, NWTSS);
//  sourceGaussian(true, AMPH, &(space2d->c[SRC5].h1), time_step, DTS, NWTSS);
//  sourceGaussian(true, AMPH, &(space2d->c[SRC5].h2), time_step, DTS, NWTSS);
#endif

#ifdef ABC_SECOND_ORDER
  abc2o2d->update();
#endif

  time_step++;
}

void CModel2D::reset()
{
  time_step = 0; // time step
  space2d->reset();
  if(abc2o2d != NULL) abc2o2d->reset();
  if(tfsf2d != NULL) tfsf2d->reset();

#ifdef SNAPSHOT
  #define SRC_X (SIZEX / 2)
  #define SRC_Y (SIZEY / 2)
  #define SRC0 (SRC_X +      SRC_Y * SIZEX)
  #define SRC1 (SRC_X + 1 +  SRC_Y * SIZEX)
  #define SRC2 (SRC_X + 2 +  SRC_Y * SIZEX)
  #define SRC3 (SRC_X + 3 +  SRC_Y * SIZEX)
  #define SRC4 (SRC_X +     (SRC_Y + 1) * SIZEX)
  #define SRC5 (SRC_X + 1 + (SRC_Y + 1) * SIZEX)
  #define SRC6 (SRC_X + 2 + (SRC_Y + 1) * SIZEX)
  #define SRC7 (SRC_X + 3 + (SRC_Y + 1) * SIZEX)
  #define SRC8 (SRC_X +     (SRC_Y + 2) * SIZEX)
  #define SRC9 (SRC_X + 1 + (SRC_Y + 2) * SIZEX)
  #define SRCA (SRC_X + 2 + (SRC_Y + 2) * SIZEX)
  #define SRCB (SRC_X + 3 + (SRC_Y + 2) * SIZEX)
  #define AMPE 0.4
  #define AMPH (AMPE / IMP0)
  space2d->c[SRC1].e = AMPE;
  space2d->c[SRC4].e = AMPE;
  space2d->c[SRC6].e = AMPE;
  space2d->c[SRC9].e = AMPE;

//  space2d->c[SRC1].h1 = AMPH;
//  space2d->c[SRC4].h1 = AMPH;
//  space2d->c[SRC5].h1 = AMPH;
//  space2d->c[SRC5].h2 = AMPH;
#endif
}

void CModel2D::inc_cut_type()
{
  if(++cut_type > 2) cut_type = 0;
}

void CModel2D::inc_field_type()
{
  if(++field_type > 4) field_type = 0;
}

void CModel2D::draw() const  // right handed space
{
  glBegin(GL_POINTS);
  #define DECADES 2.3
  #define MAX 0.5
  #define ACUT 0.03
  #define HSVH1 (-360.0 / (DECADES * log(10.0)))
  #define HSVH2 1.0 / MAX
  #define ACUTM (1.0 / (MAX * ACUT))
  // display as Ez Hx Hy
  HSV hsv;
  RGB rgb;
  size_t n;
  hsv.s = 1.0;
  hsv.v = 0.8;
  if(cut_type < 2) {
    double e, h1, h2, eh;
    for(size_t i = 0; i < SIZEX; i++) {
      for(size_t j = 0; j < SIZEY; j++) {
        n = i + j * SIZEX;

        if(field_type == 0) eh = space2d->c[n].e; // field type?
        else if(field_type == 1) eh = IMP0 * space2d->c[n].h1;
        else if(field_type == 2) eh = IMP0 * space2d->c[n].h2;
        else if(field_type == 3) {
          h1 = IMP0 * space2d->c[n].h1;
          h2 = IMP0 * space2d->c[n].h2;
          eh = sqrt(h1 * h1 + h2 * h2); // norm
        }
        else if(field_type == 4) {
          e = space2d->c[n].e;
          h1 = IMP0 * space2d->c[n].h1;
          h2 = IMP0 * space2d->c[n].h2;
          eh = sqrt(e * e + h1 * h1 + h2 * h2); // norm
        }

        double ne = fabs(eh);
        double alpha = ACUTM * ne;
        if(alpha < 0.05) continue; // skip low alphas
        hsv.h = HSVH1 * log(HSVH2 * ne);
        if(hsv.h >= 360.0) hsv.h = 359.9;
        HSVtoRGB(&hsv, &rgb);
        glColor4d(rgb.r, rgb.g, rgb.b, alpha_fac * alpha);
        double x = -0.5 + i / (SIZEX - 1.0);
        double y = -0.5 + j / (SIZEY - 1.0);
        double v = cut_type ? 0.0 : eh; // cut type surface or flat ?
        if(dither)
          if(p->run)
            glVertex3d((randpm() / SIZEX) + x,
                       (randpm() / SIZEY) + y,
                       (randpm() / SIZEZ) + v);
          else
            glVertex3d((skip * space2d->d[n] / SIZEX) + x,
                       (skip * space2d->d[2 * n] / SIZEY) + y,
                       (skip * space2d->d[3 * n] / SIZEZ) + v);
        else glVertex3d(x, y, v);
      }
    }
  }
  else {
    size_t l = (SIZEY / 2) * SIZEX;
    glColor3d(BRIGHTYELLOW);
    for(int i = 0; i < SIZEX; i++) {
      double e = space2d->c[i + l].e;
      glVertex3d(-0.5 + ((1.0 * i) / (SIZEX - 1)), 0.0, e);
    }
    glColor3d(BRIGHTCYAN);
    for(int i = 0; i < SIZEX; i++) {
      double h1 = IMP0 * space2d->c[i + l].h1;
      double h2 = IMP0 * space2d->c[i + l].h2;
      double h = sqrt(h1 * h1 + h2 * h2);
      glVertex3d(-0.5 + ((0.5 + i) / (SIZEX - 1)), h, 0.0);
    }
  }
  glEnd();

  if(display_bbox) glCallList(bbox);
}

void CModel2D::get_status(QString &s) const
{
  QString d;
  s.append("2D, ");
  if(cut_type < 2) {
    if(field_type == 0) s.append("Ez");
    else if(field_type == 1) s.append("Hx");
    else if(field_type == 2) s.append("Hy");
    else if(field_type == 3) s.append("Hxy");
    else if(field_type == 4) s.append("EzHxy");
  }
  else s.append("Ez, Hxy");
  d = ", " + QString::number(space2d->sX);
  d += "x" + QString::number(space2d->sY);
  s.append(d);
  if(dither && (cut_type != 2)) s.append(", dither");
}
