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
#include "cell3d.h"
#include "spaceEH3d.h"
#include "abc1o3d.h"
#include "tfsf3d.h"

#include "model3d.h"

#define SIZEX 61
#define SIZEY 61
#define SIZEZ 61

#define FREE_SPACE

#define PEC_SPHERE_00
//#define FIELD_PEC_SLIT

#define RICKER_PLANE
//#define GAUSSIAN_PLANE
//#define RICKER_POINT
//#define GAUSSIAN_POINT

//#define ABC_FIRST_ORDER

CModel3D::CModel3D(GLWidget *parent) : CModel(parent)
{
  space3d = NULL;
  abc1o3d = NULL;
  tfsf3d = NULL;
  objects = NULL;

  set_material(); // space & material
  reset();        // space & material
}
// ***********************************************************************
// model materials
// ***********************************************************************
void CModel3D::set_material()
{
  space3d = new CSpaceEH3d(SIZEX, SIZEY, SIZEZ);

#ifdef FREE_SPACE
  for(int i = 0; i < SIZEX * SIZEY * SIZEZ; i++) {
    space3d->c[i].cexe = space3d->c[i].ceye = space3d->c[i].ceze = 1.0; // free-space
    space3d->c[i].cexh = space3d->c[i].ceyh = space3d->c[i].cezh = DTDS3D * IMP0;
    space3d->c[i].chxh = space3d->c[i].chyh = space3d->c[i].chzh = 1.0;
    space3d->c[i].chxe = space3d->c[i].chye = space3d->c[i].chze = DTDS3D / IMP0;
  }
#endif

#ifdef PEC_SPHERE_00
  #define CR (SIZEX / 7)
  #define CX (3 * SIZEX / 4)
  #define CY (SIZEY / 2)
  #define CZ (SIZEZ / 2)
  for(int k = CZ - CR; k < CX + CR; k++) {
    for(int j = CY - CR; j < CY + CR; j++) {
      for(int i = CX - CR; i < CX + CR; i++) {
        int cr2 = CR * CR;
        if( pow((i - CX), 2) + pow((j - CY), 2) + pow((k - CZ), 2) <= cr2) {
          int n = i + j * SIZEX + k * SIZEX * SIZEY;
          space3d->c[n].cexe = 0.0;
          space3d->c[n].cexh = 0.0;
          space3d->c[n].ceye = 0.0;
          space3d->c[n].ceyh = 0.0;
          space3d->c[n].ceze = 0.0;
          space3d->c[n].cezh = 0.0;
        }
      }
    }
  }

  GLfloat mat_diffuse[] = {0.4, 0.4, 0.4, 1.0};
  GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat mat_shininess[] = {1.0};
  objects = glGenLists(1);
  float s = (float(CR) - 0.5) / SIZEX;
  float tX = 0.25;
  float tY = 0.0;
  float tZ = 0.0;
  glNewList(objects, GL_COMPILE);
    glPushMatrix();
    glTranslatef(tX, tY, tZ);
    glScalef(s, s, s);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    icosphere(5);
    glPopMatrix();
  glEndList();
#endif

  // set tfsf and abc's after material initialization!!!

#ifdef GAUSSIAN_PLANE
  #define SD 10  // tfsf aux decay size
  #define SB 3  // tfsf boundary size
  tfsf3d = new CTfsf3d(space3d, SB, SD);
#endif

#ifdef RICKER_PLANE
  #define SD 10  // tfsf aux decay size
  #define SB 3  // tfsf boundary size
  tfsf3d = new CTfsf3d(space3d, SB, SD);
#endif

#ifdef ABC_FIRST_ORDER
  abc1o3d = new CAbc1o3d(space3d);
#endif
}

// ***********************************************************************
// model field step
// ***********************************************************************
void CModel3D::step()
{
  space3d->update_h(); //  ***** update magnetic field *****

#ifdef RICKER_PLANE
  #define WTS 100
  tfsf3d->updateA();
  sourceRicker(true, 0.3 / -IMP0, tfsf3d->inpm1, time_step, WTS);
  sourceRicker(true, 0.3, tfsf3d->inp, time_step + 1, WTS);
  tfsf3d->updateB();
#endif

#ifdef GAUSSIAN_PLANE
  #define WTS 10
  #define DTS (WTS * 4) // delay time steps
  #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  tfsf3d->updateA();
  sourceGaussian(true, 0.3 / -IMP0, tfsf3d->inpm1, time_step, DTS, NWTSS);
  sourceGaussian(true, 0.3, tfsf3d->inp, time_step + 1, DTS, NWTSS);
  tfsf3d->updateB();
#endif

  space3d->update_e();  // ***** update electric field *****

#ifdef RICKER_POINT
  #define WTS 200
  #define SRC_X (SIZEX / 2)
  #define SRC_Y (SIZEY / 2)
  #define SRC_Z (SIZEZ / 2)
  #define SRC (SRC_X + SRC_Y * SIZEX + SRC_Z * SIZEX * SIZEY)
  sourceRicker(false, 10.0, &(space3d->c[SRC].ez), time_step, WTS);
#endif

#ifdef GAUSSIAN_POINT
  #define WTS 10
  #define DTS (WTS * 4) // delay time steps
  #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  #define SRC_X (SIZEX / 4)
  #define SRC_Y (SIZEY / 2)
  #define SRC_Z (SIZEZ / 2)
  #define SRC (SRC_X + SRC_Y * SIZEX + SRC_Z * SIZEX * SIZEY)
  sourceGaussian(false, 10.0, &(space3d->c[SRC].ez), time_step, DTS, NWTSS);
#endif

#ifdef ABC_FIRST_ORDER
  abc1o3d->update_e();
#endif

#ifdef FIELD_PEC_SLIT
  for(int k = 0; k < 50; k++) {
    for(int j = 0; j < SIZEY; j++) {
      int n = 70 + j * SIZEX + k * SIZEX * SIZEY;
      space3d->c[n].ex = space3d->c[n].ey = space3d->c[n].ez = 0.0;
    }
  }
  for(int k = 51; k < SIZEZ; k++) {
    for(int j = 0; j < SIZEY; j++) {
        int n = 70 + j * SIZEX + k * SIZEX * SIZEY;
        space3d->c[n].ex = space3d->c[n].ey = space3d->c[n].ez = 0.0;
    }
  }
#endif

  time_step++;
}

void CModel3D::reset()
{
  time_step = 0; // time step
  alpha_fac = 1.0;
  space3d->reset();
  if(abc1o3d != NULL) abc1o3d->reset();
  if(tfsf3d != NULL) tfsf3d->reset();
}

void CModel3D::inc_cut_type()
{
  if(++cut_type > 4) cut_type = 0;
}

void CModel3D::inc_field_type()
{
  if(++field_type > 8) field_type = 0;
}

void CModel3D::draw() const  // right handed space
{
  if(objects != 0) glCallList(objects);

  glBegin(GL_POINTS);
  #define DECADES 2.3
  #define MAX 0.5
  #define HSVH1 (-360.0 / (DECADES * log(10.0)))
  #define HSVH2 1.0 / MAX
  #define ACUT 0.03 // alpha-one cut-off factor
  #define ACUTM (1.0 / (MAX * ACUT)) // alpha cut multiplier
  HSV hsv;
  RGB rgb;
  long ii, jj, kk, maxii, maxjj, maxkk;
  size_t n;
  double eh, ex, ey, ez, hx, hy, hz;
  double x, y, z;
  hsv.s = 1.0;
  hsv.v = 0.8;

  if(face == 3) {
    maxkk = SIZEZ;
    maxjj = SIZEY;
    maxii = SIZEX;
  }
  else if(face == 2) {
    maxkk = SIZEY;
    maxjj = SIZEX;
    maxii = SIZEZ;
  }
  else{
    maxkk = SIZEX;
    maxjj = SIZEZ;
    maxii = SIZEY;
  }

  long i = 0;
  long j = 0;
  long k = 0;
  long sxmax, symax, szmax, smin;
  sxmax = symax = szmax = smin = 0;
  if(!display_boundary && (tfsf3d != NULL)) {
    sxmax = SIZEX - tfsf3d->sb;
    symax = SIZEY - tfsf3d->sb;
    szmax = SIZEZ - tfsf3d->sb;
    smin = tfsf3d->sb;
  }
  for(kk = 0; kk < maxkk; kk += skip) {
    if(face == 3)          k = zp ? kk : maxkk - 1 - kk;
    else if(face == 2)     j = yp ? kk : maxkk - 1 - kk;
    else                   i = xp ? kk : maxkk - 1 - kk;
    for(jj = 0; jj < maxjj; jj += skip) {
      if(face == 3)          j = yp ? jj : maxjj - 1 - jj;
      else if(face == 2)     i = xp ? jj : maxjj - 1 - jj;
      else                   k = zp ? jj : maxjj - 1 - jj;
      for(ii = 0; ii < maxii; ii += skip) {
        if(face == 3)          i = xp ? ii : maxii - 1 - ii;
        else if(face == 2)     k = zp ? ii : maxii - 1 - ii;
        else                   j = yp ? ii : maxii - 1 - ii;

        if((cut_type == 1) && j < (SIZEY / 2)) continue; // display half: skip
        if((cut_type == 2) && j != (SIZEY / 2)) continue; // display slice: skip
        if((cut_type == 3) && ( k != SIZEZ / 2)) continue; // display surface: skip
        if((cut_type == 4) && (j != (SIZEY / 2) || k != SIZEZ / 2)) continue; // display line: skip
        if(sxmax != 0) if(i >= sxmax) continue;
        if(symax != 0) if(j >= symax) continue;
        if(szmax != 0) if(k >= szmax) continue;
        if(smin != 0) if((i < smin) || (j < smin) || (k < smin)) continue;

        n = i + j * SIZEX + k * SIZEX * SIZEY;
        eh = 0.0; // play it safe

        if(field_type == 1) eh = space3d->c[n].ex; // field type?
        else if(field_type == 2) eh = space3d->c[n].ey;
        else if(field_type == 3) eh = space3d->c[n].ez;
        else if(field_type == 0) {
          ex = space3d->c[n].ex;
          ey = space3d->c[n].ey;
          ez = space3d->c[n].ez;
          eh = sqrt(ex * ex + ey * ey + ez * ez); // norm
        }
        else if(field_type == 5) eh = IMP0 * space3d->c[n].hx;
        else if(field_type == 6) eh = IMP0 * space3d->c[n].hy;
        else if(field_type == 7) eh = IMP0 * space3d->c[n].hz;
        else if(field_type == 5) {
          hx = IMP0 * space3d->c[n].hx;
          hy = IMP0 * space3d->c[n].hy;
          hz = IMP0 * space3d->c[n].hz;
          eh = sqrt(hx * hx + hy * hy + hz * hz); // norm
        }
        else if(field_type == 8) {
          ex = space3d->c[n].hx;
          ey = space3d->c[n].hy;
          ez = space3d->c[n].hz;
          hx = IMP0 * space3d->c[n].hx;
          hy = IMP0 * space3d->c[n].hy;
          hz = IMP0 * space3d->c[n].hz;
          eh = sqrt(ex * ex + ey * ey + ez * ez + hx * hx + hy * hy + hz * hz); // norm
        }

        double alpha = ACUTM * fabs(eh);
        if((cut_type < 3) && (alpha < 0.05)) continue; // skip low alphas
        if(cut_type == 4) glColor4d(BLACK, 1.0); // line?
        else {
          hsv.h = HSVH1 * log(HSVH2 * fabs(eh)); // set color
          if(hsv.h >= 360.0) hsv.h = 359.9; // clamp at max
          HSVtoRGB(&hsv, &rgb);
          glColor4d(rgb.r, rgb.g, rgb.b, alpha_fac * alpha);
        }
        x = -0.5 + i / (SIZEX - 1.0);
        y = -0.5 + j / (SIZEY - 1.0);
        if(cut_type > 2) z = eh; // surface or line?
        else z = -0.5 + k / (SIZEZ - 1.0);
        if(dither && (cut_type != 4))
          if(p->run)
            glVertex3d((skip * randpm() / SIZEX) + x,
                       (skip * randpm() / SIZEY) + y,
                       (skip * randpm() / SIZEZ) + z);
          else
            glVertex3d((skip * space3d->d[n] / SIZEX) + x,
                       (skip * space3d->d[2 * n] / SIZEY) + y,
                       (skip * space3d->d[3 * n] / SIZEZ) + z);
        else glVertex3d(x, y, z);
      }
    }
  }
  glEnd();

  if(display_bbox) glCallList(bbox);

  if(display_yee) {
    glPushMatrix();
    glScaled(1.0 / SIZEX, 1.0 / SIZEY, 1.0 / SIZEZ);
    for(size_t k = 0; k <= SIZEZ; k++) {
      for(size_t j = 0; j <= SIZEY; j++) {
        for(size_t i = 0; i <= SIZEX; i++) {
          glPushMatrix();
          glTranslated(i - SIZEX / 2.0, j - SIZEY / 2.0, k - SIZEZ / 2.0);
          space3d->draw_yee();
          glPopMatrix();
        }
      }
    }
    glPopMatrix();
  }
}

void CModel3D::get_status(QString &s) const
{
  QString d;
  s.append("3D, ");
  if(field_type == 0) s.append("Exyz");
  else if(field_type == 1) s.append("Ex");
  else if(field_type == 2) s.append("Ey");
  else if(field_type == 3) s.append("Ez");
  else if(field_type == 4) s.append("Hxyz");
  else if(field_type == 5) s.append("Hx");
  else if(field_type == 6) s.append("Hy");
  else if(field_type == 7) s.append("Hz");
  else if(field_type == 8) s.append("ExyzHxyz");
  d = ", " + QString::number(space3d->sX);
  d += "x" + QString::number(space3d->sY);
  d += "x" + QString::number(space3d->sZ);
  s.append(d);
  if(cut_type == 0) s.append(", full");
  else if(cut_type == 1) s.append(", half");
  else if(cut_type == 2) s.append(", slice");
  else if(cut_type == 3) s.append(", surface");
  else if(cut_type == 4) s.append(", line");
  if(dither && (cut_type != 4)) s.append(", dither");
}
