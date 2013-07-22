/*
GL_10
An OpenGL+Qt4 FDTD electromagnetic simulation & visualization program.

Copyright (C) 2005-2013 John Rugis

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
#include <math.h>

#include "defs.h"
#include "utils.h"
#include "source.h"
#include "cell1d.h"
#include "spaceEH1d.h"
#include "abc1o1d.h"
#include "abc2o1d.h"

#include "model1d.h"

// NOTES:
//  1) LHS e-field and RHS h-field not updated

#define SIZEX 100

#define FREE_SPACE
//#define LOSSY_E_SPACE
//#define LOSSLESS_DIELECTRIC_SPACE
//#define HALF_SPACE_LOSSY_E
//#define HALF_SPACE_LOSSLESS_DIELECTRIC
//#define HALF_SPACE_LOSSY_DIELECTRIC
//#define HALF_SPACE_LOSSLESS_DIELECTRIC_MATCHED_LOSSY_RHS

//#define GAUSSIAN_LHS
#define GAUSSIAN_TFSF
//#define GAUSSIAN_INTERNAL
//#define IMPULSE_LHS
//#define IMPULSE_INTERNAL
//#define IMPULSE_TFSF
//#define IMPULSE_TO_RIGHT
//#define IMPULSE_OUTWARD
//#define SINE_LHS
//#define SINE_INTERNAL
//#define SINE_TFSF
//#define RICKER_LHS
//#define RICKER_TFSF

//#define FIELD_PEC
//#define FIELD_PMC

//#define ABC_SIMPLE_LHS
//#define ABC_SIMPLE_RHS
//#define ABC_FIRST_ORDER
//#define ABC_SECOND_ORDER

CModel1D::CModel1D(GLWidget *parent) : CModel(parent)
{
  space1d = NULL; // space & material
  abc1o1d = NULL; // advanced abc's
  abc2o1d = NULL;

  set_material(); // space & material
  reset();        // space & material
}

// ***********************************************************************
// model materials (initial)
// ***********************************************************************
void CModel1D::set_material()
{
  space1d = new CSpaceEH1d(SIZEX);

#ifdef FREE_SPACE
  for(int i = 0; i < SIZEX; i++) {
    space1d->c[i].cee = 1.0;  // free-space
    space1d->c[i].ceh = IMP0;
    space1d->c[i].chh = 1.0;
    space1d->c[i].che = 1.0 / IMP0;
  }
#endif

#ifdef LOSSY_E_SPACE
  #define LOSS 0.01
  for(int i = 0; i < SIZEX; i++) {
    space1d->c[i].cee = (1.0 - LOSS) / (1.0 + LOSS); // e-field lossy
    space1d->c[i].ceh = IMP0 / (1.0 + LOSS);
    space1d->c[i].chh = 1.0;
    space1d->c[i].che = 1.0 / IMP0;
  }
#endif

#if defined LOSSLESS_DIELECTRIC_SPACE
  #define EPSR 2.0
  for(int i = 0; i <  SIZEX; i++) {
    space1d->c[i].cee = 1.0;
    space1d->c[i].ceh = IMP0 / EPSR;
    space1d->c[i].chh = 1.0;  // free-space
    space1d->c[i].che = 1.0 / IMP0;
  }
#endif

#ifdef HALF_SPACE_LOSSY_E
  #define HS (SIZEX / 2)
  #define HS_LOSS 0.03
  for(int i = 0; i < SIZEX; i++) {
    if(i < HS) {              // free-space
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0;
    }
    else {                   // e-field lossy
      space1d->c[i].cee = (1.0 - HS_LOSS) / (1.0 + HS_LOSS);
      space1d->c[i].ceh = IMP0 / (1.0 + HS_LOSS);
    }
  }
  for(int i = 0; i < SIZEX; i++) {
    space1d->c[i].chh = 1.0;  // free-space
    space1d->c[i].che = 1.0 / IMP0;
  }
#endif

#ifdef HALF_SPACE_LOSSLESS_DIELECTRIC
  #define HS (SIZEX / 2)
  #define RHS_EPSR 2.0  // right half-space
  for(int i = 0; i <  SIZEX; i++) {
    if(i<HS) {              // free-space
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0;
    }
    else {   // lossless dielectric
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0 / RHS_EPSR;
    }
  }
  for(int i = 0; i < SIZEX; i++) {
    space1d->c[i].chh = 1.0;  // free-space
    space1d->c[i].che = 1.0 / IMP0;
  }
#endif

#if defined HALF_SPACE_LOSSY_DIELECTRIC
  #define HS (0.5 * SIZEX)
  #define RHS_LOSS 0.03
  #define RHS_EPSR 2.0
  for(int i = 0; i < SIZEX; i++) {
    if(i<HS) {              // free-space
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0;
    }
    else {   // lossy dielectric
      space1d->c[i].cee = (1.0 - RHS_LOSS) / (1.0 + RHS_LOSS);
      space1d->c[i].ceh = IMP0 / RHS_EPSR / (1.0 + RHS_LOSS);
    }
  }
  for(int i = 0; i < SIZEX; i++) {
    space1d->c[i].chh = 1.0;  // free-space
    space1d->c[i].che = 1.0 / IMP0;
  }
#endif

#ifdef HALF_SPACE_LOSSLESS_DIELECTRIC_MATCHED_LOSSY_RHS
  #define HS (SIZEX / 2)
  #define RBL_START (4 * (SIZEX / 5))
  #define RHS_EPSR 9.0
  #define RBL_LOSS 0.01
  for(int i = 0; i < SIZEX; i++) {
    if(i < HS) {              // free-space
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0;
    }
    else if(i < RBL_START) {   // lossless dielectric
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0 / RHS_EPSR;
    }
    else {                   // e-field matched lossy boundary layer
      space1d->c[i].cee = (1.0 - RBL_LOSS) / (1.0 + RBL_LOSS);
      space1d->c[i].ceh = IMP0 / RHS_EPSR / (1.0 + RBL_LOSS);
    }
  }
  for(int i=0; i < SIZEX; i++) {
    if(i < RBL_START) {     // free-space
      space1d->c[i].chh = 1.0;
      space1d->c[i].che = 1.0 / IMP0;
    }
    else {                // h-field matched lossy boundary layer
      space1d->c[i].chh = (1.0 - RBL_LOSS) / (1.0 + RBL_LOSS);
      space1d->c[i].che = 1.0 / IMP0 / (1.0 + RBL_LOSS);
    }
  }
#endif

// set abc's after material initialization!!!

#ifdef ABC_FIRST_ORDER
  abc1o1d = new CAbc1o1d(space1d);
#endif
#ifdef ABC_SECOND_ORDER
  abc2o1d = new CAbc2o1d(space1d);
#endif

}

// ***********************************************************************
// model field step
// ***********************************************************************
void CModel1D::step()
{
  space1d->update_h(); // ***** update magnetic field *****

#ifdef GAUSSIAN_LHS
  #define WTS 30
  #define DTS (WTS * 4) // delay time steps
  #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  #ifdef D22
    sourceGaussian(false, 0.5, &(space1d->c[0].e), time_step, DTS, NWTSS);
  #endif
  #ifdef D24
    sourceGaussian(false, 0.5, &(space1d->c[0].e), time_step + 1, DTS, NWTSS);
    sourceGaussian(false, 0.5, &(space1d->c[1].e), time_step, DTS, NWTSS);
  #endif
#endif

#ifdef GAUSSIAN_TFSF
  #define SRC_POS (SIZEX / 6)   // source position
  #define WTS 30
  #define DTS (WTS * 4) // delay time steps
  #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  sourceGaussian(true,  0.4 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, DTS, NWTSS); // free-space TFSF source
  sourceGaussian(true, 0.4, &(space1d->c[SRC_POS].e), time_step + 1, DTS, NWTSS);
#endif

#ifdef GAUSSIAN_INTERNAL
  #define SRC_POS (SIZEX / 3)   // source position
  #define WTS 30
  #define DTS (WTS * 4) // delay time steps
  #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  sourceGaussian(true, 0.4, &(space1d->c[SRC_POS].e), time_step, DTS, NWTSS);
  sourceGaussian(true, 0.4, &(space1d->c[SRC_POS + 1].e), time_step, DTS, NWTSS);
#endif

#ifdef SINE_LHS
  #define TSW (SIZEX / 4.0) // time steps per wavelength
  #define OMEGA (2.0 * PI / TSW) // omega
  sourceSine(false, 0.25, &(space1d->c[0].e), time_step, OMEGA);
#endif

#ifdef SINE_INTERNAL
  #define SRC_POS (SIZEX / 6)   // source position
  #define TSW (SIZEX / 4.0) // time steps per wavelength
  #define OMEGA (2.0 * PI / TSW) // omega
  sourceSine(true, 0.4, &(space1d->c[SRC_POS].e), time_step, OMEGA);
  sourceSine(true, 0.4, &(space1d->c[SRC_POS + 1].e), time_step, OMEGA);
#endif

#ifdef SINE_TFSF
  #define SRC_POS (SIZEX / 6)   // source position
  #define TSW (SIZEX / 5.0) // time steps per wavelength
  #define OMEGA (2.0 * PI / TSW) // omega
  sourceSine(true, 0.25 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, OMEGA);  // free-space TFSF source
  sourceSine(true, 0.25, &(space1d->c[SRC_POS].e), time_step + 1, OMEGA);
#endif

#ifdef IMPULSE_LHS
  #define ON 5
  #define OFF 6
  sourceImpulse(false, 0.4, &(space1d->c[0].e), time_step, ON, OFF);
#endif

#ifdef IMPULSE_INTERNAL
  #define SRC_POS (SIZEX / 3)   // source position
  #define ON 5
  #define OFF 10
  sourceImpulse(true, 0.4, &(space1d->c[SRC_POS].e), time_step, ON, OFF);
  sourceImpulse(true, 0.4, &(space1d->c[SRC_POS + 1].e), time_step, ON, OFF);
#endif

#ifdef IMPULSE_TFSF
  #define SRC_POS (SIZEX / 6)   // source position
  #define ON 5
  #define OFF 45
  sourceImpulse(true, 0.4 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, ON, OFF); // free-space TFSF source
  sourceImpulse(true, 0.4, &(space1d->c[SRC_POS].e), time_step + 1, ON, OFF);
#endif

#ifdef RICKER_LHS
  #define WTS 200
  sourceRicker(false, 0.4, &(space1d->c[0].e), time_step, WTS);
#endif

#ifdef RICKER_TFSF
  #define SRC_POS (SIZEX / 6)   // source position
  #define WTS 200
  #define TSW (SIZEX / 10.0) // time steps per wavelength
  #define OMEGA (2.0 * PI / TSW) // omega
  sourceRicker(true, 0.4 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, WTS);  // free-space TFSF source
  sourceRicker(true, 0.4, &(space1d->c[SRC_POS].e), time_step + 1, WTS);
#endif

#ifdef FIELD_PMC
  #define POS (5 * (SIZEX / 6))
  space1d->c[POS].h = 0.0; // embedded electricmagnetic conductor
#endif
#ifdef ABC_SIMPLE_LHS
  space1d->c[0].e = space1d->c[1].e; // LHS - simple ABC
#endif
#ifdef ABC_FIRST_ORDER
  abc1o1d->update_h();
#endif
#ifdef ABC_SECOND_ORDER
  abc2o1d->update_h();
#endif

  space1d->update_e();  // ***** update electric field *****

#ifdef ABC_SIMPLE_RHS
  space1d->c[SIZEX - 1].h = space1d->c[SIZEX - 2].h; // simple RHS h-field ABC
#endif
#ifdef ABC_FIRST_ORDER
  abc1o1d->update_e();
#endif
#ifdef ABC_SECOND_ORDER
  abc2o1d->update_e();
#endif
#ifdef FIELD_PEC
  #define POS (5 * (SIZEX / 6))
  space1d->c[POS].e = 0.0; // embedded electric conductor
#endif

  time_step++;
}

void CModel1D::reset()
{
  time_step = 0; // time step
  space1d->reset();
  if(abc1o1d != NULL) abc1o1d->reset();
  if(abc2o1d != NULL) abc2o1d->reset();

#ifdef IMPULSE_TO_RIGHT
  #define SRC_POS (SIZEX / 3)
  #define AMP 0.4
  #define WIDTH 10
  for(int i = 0; i < WIDTH; i++) {
    space1d->c[SRC_POS + i].e = AMP;
    space1d->c[SRC_POS + i - 1].h = -AMP / IMP0;
  }
#endif

#ifdef IMPULSE_OUTWARD
  #define SRC_POS (SIZEX / 3)
  #define AMP 0.4
  space1d->c[SRC_POS].e = AMP;
  space1d->c[SRC_POS + 1].e = AMP;
#endif

}

void CModel1D::inc_cut_type()
{
}

void CModel1D::inc_field_type()
{
  if(++field_type > 3) field_type = 0;
}

void CModel1D::draw() const  // right handed space
{
  glBegin(GL_POINTS);
  if(field_type != 2) { // display Ez
    glColor3d(BRIGHTYELLOW);
    for(int i = 0; i < SIZEX; i++)
      glVertex3d(-0.5 + (1.0 * i / SIZEX), 0.0, space1d->c[i].e);
    if(field_type == 1) {
      glColor3d(BRIGHTRED);
      for(int i = 0; i < SIZEX; i++) {
          glVertex3d(-0.5 + (1.0 * i / SIZEX), 0.0, space1d->eMax[i]);
          glVertex3d(-0.5 + (1.0 * i / SIZEX), 0.0, space1d->eMin[i]);
      }
    }
  }
  if(field_type >= 2) { // display Hy
    glColor3d(BRIGHTCYAN);
    for(int i = 0; i < SIZEX; i++)
      glVertex3d(-0.5 + ((0.5 + i) / SIZEX), IMP0 * space1d->c[i].h, 0.0);
  }
  glEnd();

  if(display_bbox) glCallList(bbox);
}

void CModel1D::get_status(QString &s) const
{
  QString d;
  s.append("1D, ");
  if(field_type == 0) s.append("Ez");
  else if(field_type == 1) s.append("Ez-peak");
  else if(field_type == 2) s.append("Hy");
  else if(field_type == 3) s.append("EzHy");
  d = ", " + QString::number(space1d->size);
  s.append(d);
}


