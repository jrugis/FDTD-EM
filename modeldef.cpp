/*

GL_10
An OpenGL+Qt4 FDTD electromagnetic simulation & visualization program.

portions Copyright (C) 2005-2011 J.Rugis
portions Copyright (C) 2012 Institute of Earth Science & Engineering

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

j.rugis@auckland.ac.nz
*/

#include <fstream>
#include <math.h>

#include "defs.h"
#include "utils.h"
#include "glwidget.h"
#include "source.h"
#include "cell1d.h"
#include "cell2d.h"
#include "cell3d.h"
#include "spaceEH1d.h"
#include "spaceEH2d.h"
#include "spaceEH3d.h"
#include "abc1o1d.h"
#include "abc1o3d.h"
#include "abc2o1d.h"
#include "abc2o2d.h"
#include "tfsf2d.h"
#include "tfsf3d.h"

#include "model.h"

#define M1000
#define F1000

#if defined M0000 || defined M0001 || defined M0002 || defined M0003 || defined M0004
  #define SIZEX 300
  #define DECADES 2.3
  #define MAX 0.5
  #define ACUT 0.03
  #define HSVH1 (-360.0 / (DECADES * log(10.0)))
  #define HSVH2 1.0 / MAX
  #define ACUTM (1.0 / (MAX * ACUT))
#endif

#if defined M1000
  #define SIZEX 401
  #define SIZEY 401
  #define SIZEZ 401
#endif

#if defined M2000
  #define SIZEX 101
  #define SIZEY 101
  #define SIZEZ 101
  #define DECADES 5.0
  #define MAX 0.2
  #define HSVH1 (-360.0 / (DECADES * log(10.0)))
  #define HSVH2 1.0 / MAX
  #define ACUT 0.0003 // alpha-one cut-off factor
  #define ACUTM (1.0 / (MAX * ACUT)) // alpha cut multiplier
#endif

// ***********************************************************************
// model materials
// ***********************************************************************
void CModel::set_material()
{
#if defined M0000
  space1d = new CSpaceEH1d(SIZEX);
  for(int i = 0; i < SIZEX; i++) {
    space1d->c[i].cee = 1.0;  // free-space
    space1d->c[i].ceh = IMP0;
    space1d->c[i].chh = 1.0;
    space1d->c[i].che = 1.0 / IMP0;
  }
  // set abc's after material initialization!!!
  //abc1o1d = new CAbc1o1d(space1d);
  abc2o1d = new CAbc2o1d(space1d);
#endif
#if defined M0001
 #define HS (0.5 * SIZEX)
 #define RHS_EPSR 9.0  // right half-space
  space1d = new CSpaceEH1d(SIZEX);
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
  // set abc's after material initialization!!!
  //abc1o1d = new CAbc1o1d(space1d);
  abc2o1d = new CAbc2o1d(space1d);
#endif
#if defined M0002
 #define HS (0.5 * SIZEX)
 #define RHS_LOSS 0.03
 #define RHS_EPSR 4.0
  space1d = new CSpaceEH1d(SIZEX);
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
  // set abc's after material initialization!!!
  //abc1o1d = new CAbc1o1d(space1d);
  abc2o1d = new CAbc2o1d(space1d);
#endif
#if defined M0003
 #define HS (0.5 * SIZEX)
 #define RBL_START (0.9 * SIZEX)
 #define RHS_EPSR 9.0
 #define RBL_LOSS 0.02
  space1d = new CSpaceEH1d(SIZEX);
  for(int i = 0; i < SIZEX; i++) {
    if(i < HS) {              // free-space
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0;
    }
    else if(i < RBL_START) {   // lossless dielectric
      space1d->c[i].cee = 1.0;
      space1d->c[i].ceh = IMP0 / RHS_EPSR;
    }
    else {                   // matched lossy boundary layer
      space1d->c[i].cee = (1.0 - RBL_LOSS) / (1.0 + RBL_LOSS);
      space1d->c[i].ceh = IMP0 / RHS_EPSR / (1.0 + RBL_LOSS);
    }
  }
  for(int i=0; i < SIZEX; i++) {
    if(i < RBL_START) {     // free-space
      space1d->c[i].chh = 1.0;
      space1d->c[i].che = 1.0 / IMP0;
    }
    else {                // matched lossy boundary layer
      space1d->c[i].chh = (1.0 - RBL_LOSS) / (1.0 + RBL_LOSS);
      space1d->c[i].che = 1.0 / IMP0 / (1.0 + RBL_LOSS);
    }
  }
#endif
#if defined M0004
 #define HS (0.5 * SIZEX)
 #define RHS_EPSR 9.0
  space1d = new CSpaceEM1d(SIZEX);
  for(int i = 0; i < SIZEX; i++) {
    if(i < HS) {              // free-space
      space1d->ceze[i] = 1.0;
      space1d->cezh[i] = IMP0;
    }
    else {   // lossless dielectric
      space1d->ceze[i] = 1.0;
      space1d->cezh[i] = IMP0 / RHS_EPSR;
    }
  }
  for(int i = 0; i < SIZEX; i++) {
    space1d->chyh[i] = 1.0;  // free-space
    space1d->chye[i] = 1.0 / IMP0;
  }
  // set abc's after material initialization!!!
  //abc1o1d = new CAbc1o1d(space1d);
  abc2o1d = new CAbc2o1d(space1d);
#endif
#if defined M1000
  #define RHS_EPSR 8.0
  space2d = new CSpaceEH2d(SIZEX, SIZEY);
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
  /*// PEC's
  #define CR 50
  #define CX 100
  #define CY 150
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
  #define CR 40
  #define CX 200
  #define CY 200
  for(int j = CY - CR; j < CY + CR; j++) {
    for(int i = CX - CR; i < CX + CR; i++) {
      int cr2 = CR * CR;
      if( pow((i - CX), 2) + pow((j - CY), 2) <= cr2) {
        int n = i + j * SIZEX;
        space2d->c[n].cee = 0.0;
        space2d->c[n].ceh = 0.0;
      }
    }
  }*/
  for(int j = 0; j < 100; j++) {
    space2d->c[250 + (j +150) * SIZEX].cee = 0.0;
    space2d->c[250 + (j +150) * SIZEX].ceh = 0.0;
    space2d->c[350 + (j +150) * SIZEX].cee = 0.0;
    space2d->c[350 + (j +150) * SIZEX].ceh = 0.0;
    space2d->c[(j + 250) + 150 * SIZEX].cee = 0.0;
    space2d->c[(j + 250) + 150 * SIZEX].ceh = 0.0;
    space2d->c[(j + 250) + 250 * SIZEX].cee = 0.0;
    space2d->c[(j + 250) + 250 * SIZEX].ceh = 0.0;
  }
  //for(int j = 100; j < 200; j++) {
  //  space2d->c[200 + j * SIZEX].cee = 0.0;
  //  space2d->c[200 + j * SIZEX].ceh = 0.0;
  //}
  // set tfsf and abc's after material initialization!!!
  #define SB 20   // tfsf boundary size
  #define SD 20  // aux boundary size
  tfsf2d = new CTfsf2d(space2d, SB, SD);
  abc2o2d = new CAbc2o2d(space2d);
#endif
#if defined M2000
  space3d = new CSpaceEH3d(SIZEX, SIZEY, SIZEZ);
  for(int i = 0; i < SIZEX * SIZEY * SIZEZ; i++) {
    space3d->c[i].cexe = space3d->c[i].ceye = space3d->c[i].ceze = 1.0; // free-space
    space3d->c[i].cexh = space3d->c[i].ceyh = space3d->c[i].cezh = DTDS3D * IMP0;
    space3d->c[i].chxh = space3d->c[i].chyh = space3d->c[i].chzh = 1.0;
    space3d->c[i].chxe = space3d->c[i].chye = space3d->c[i].chze = DTDS3D / IMP0;
  }
  // PEC's
  #define CR 18
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
  // set tfsf and abc's after material initialization!!!
  #define SB 3  // tfsf boundary size
  #define SD 10  // tfsf aux decay size
  tfsf3d = new CTfsf3d(space3d, SB, SD);
  abc1o3d = new CAbc1o3d(space3d);
#endif
}

// ***********************************************************************
// model field step
// ***********************************************************************
void CModel::step()
{
#if defined F0000
    #define WTS 30
    #define DTS (WTS * 4) // delay time steps
    #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
  space1d->update_h(); // update magnetic field
  sourceGaussian(false, 0.4, &(space1d->c[0].e), time_step, DTS, NWTSS);
  space1d->update_e();  // update electric field
  space1d->c[SIZEX - 1].h = space1d->c[SIZEX - 2].h; // simple RHS h-field ABC
#endif
#if defined F0001
    #define WTS 30
    #define DTS (WTS * 4) // delay time steps
    #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
    #define SRC_POS (SIZEX / 6)   // source position
  space1d->update_h(); // update magnetic field
  sourceGaussian(true,  0.4 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, DTS, NWTSS); // free space TFSF source
  sourceGaussian(true, 0.4, &(space1d->c[SRC_POS].e), time_step + 1, DTS, NWTSS);
  space1d->update_e();  // update electric field
#endif
#if defined F0002
    #define TSW (SIZEX / 10.0) // time steps per wavelength
    #define OMEGA (2.0 * PI / TSW) // omega
    //#define WTS 30
    //#define DTS (WTS * 4) // delay time steps
    //#define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
    #define SRC_POS (SIZEX / 6)   // source position
  space1d->update_h(); // update magnetic field
  //sourceGaussian(true,  0.4 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, DTS, NWTSS); // free space TFSF source
  //sourceGaussian(true, 0.4, &(space1d->c[SRC_POS].e), time_step + 1, DTS, NWTSS);
  sourceSine(true, 0.25 / -IMP0, &(space1d->c[SRC_POS - 1].h), time_step, OMEGA); // TFSF source
  sourceSine(true, 0.25, &(space1d->c[SRC_POS].e), time_step + 1, OMEGA); // additive source
  space1d->c[0].e = space1d->c[1].e; // LHS - simple ABC
  space1d->update_e();  // update electric field
  space1d->c[SIZEX - 1].h = space1d->c[SIZEX - 2].h; // RHS - simple ABC
#endif
#if defined F0003
    #define WTS 30
    //#define WTS 200
    #define DTS (WTS * 4) // delay time steps
    #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
    //#define SRC_POS (0.05 * SIZEX)   // source position
    #define SRC_POS (SIZEX / 3)   // source position
  space1d->update_h(); // ****** update magnetic field ******
  //space1d->hy[200] = 0.0; // embedded magnetic conductor
  //sourceGaussian(true,  0.4 / -IMP0, space1d->hy, time_step, SRC_POS - 1, DTS, NWTSS); // free space TFSF source
  //sourceGaussian(true, 0.4, space1d->ez, time_step + 1, SRC_POS, DTS, NWTSS);
  sourceGaussian(true, 0.4, &(space1d->c[SRC_POS].e), time_step, DTS, NWTSS);
  //sourceRicker(true, 0.2 / -IMP0, space1d->hy, time_step, SRC_POS - 1, WTS);
  //sourceRicker(true, 0.2, space1d->ez, time_step + 1, SRC_POS, WTS);
  //sourceRicker(true, 0.4, space1d->ez, time_step, SRC_POS, WTS);
  //abc1o1d->update_h();
  //abc2o1d->update_h();
  space1d->update_e();  // ****** update electric field ******
  //space1d->ez[200] = 0.0; // embedded electric conductor
  //abc1o1d->update_e();
  //abc2o1d->update_e();
#endif
#if defined F0004
    #define TSW (SIZEX / 8.0) // time steps per wavelength
    #define OMEGA (2.0 * PI / TSW) // omega
    #define SRC_POS (0.15 * SIZEX)   // source position
  space1d->update_h(); // ****** update magnetic field ******
  sourceSine(true, 0.25 / -IMP0, space1d->hy, time_step, SRC_POS - 1, OMEGA); // TFSF source
  sourceSine(true, 0.25, space1d->ez, time_step + 1, SRC_POS, OMEGA); // additive source
  //foabc->update_h();
  abc2o1d->update_h();
  space1d->update_e();  // ****** update electric field ******
  //foabc->update_e();
  abc2o1d->update_e();
#endif
#if defined F1000
    #define WTS 200
    #define DTS (WTS * 4) // delay time steps
    #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
    #define SRC_X (SIZEX / 4)
    #define SRC_Y (SIZEY / 2)
    #define SRC (SRC_X + SRC_Y * SIZEX)
  space2d->update_h(); // update magnetic field
  tfsf2d->updateA();
  sourceRicker(true, 0.3 / -IMP0, tfsf2d->inpm1, time_step, WTS);
  sourceRicker(true, 0.3, tfsf2d->inp, time_step + 1, WTS);
  tfsf2d->updateB();
  space2d->update_e();  // update electric field
  //sourceRicker(false, 0.8, &(space2d->c[SRC].e), time_step, WTS);
  //sourceGaussian(false, 1.0, &(space2d->c[SRC].e), time_step, DTS, NWTSS);
  abc2o2d->update();
#endif
#if defined F2000
    #define WTS 100
    #define DTS (WTS * 4) // delay time steps
    #define NWTSS (-1.0 * pow(WTS, 2)) // negative ((width time steps) squared)
    #define SRC_X (SIZEX / 2)
    #define SRC_Y (SIZEY / 2)
    #define SRC_Z (SIZEZ / 2)
    #define SRC (SRC_X + SRC_Y * SIZEX + SRC_Z * SIZEX * SIZEY)
  space3d->update_h(); // update magnetic field
  tfsf3d->updateA();
  sourceRicker(true, 0.3 / -IMP0, tfsf3d->inpm1, time_step, WTS);
  sourceRicker(true, 0.3, tfsf3d->inp, time_step + 1, WTS);
  tfsf3d->updateB();
  space3d->update_e();  // update electric field
  //sourceRicker(false, 1.0, &(space3d->c[SRC].ez), time_step, WTS);
  abc1o3d->update_e();
  // ***** PEC slit
  /*for(int k = 0; k < 50; k++) {
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
  }*/
#endif
  time_step++;
}
