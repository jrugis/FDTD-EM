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

#include <QMessageBox>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>

#include "utils.h"

void fatalError(QString message)
{
  QMessageBox::information(0, "Fatal Error", message);
  exit(-1);
}

void normalize(float v[3])
{
  GLfloat d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d == 0.0) {
    fatalError("Normalize zero length vector.");
    return;
  }
  v[0] /= d; v[1] /= d; v[2] /= d;
}

void drawtriangle(float *v1, float *v2, float *v3)
{
  glBegin(GL_TRIANGLES);
    glNormal3fv(v1); glVertex3fv(v1);
    glNormal3fv(v2); glVertex3fv(v2);
    glNormal3fv(v3); glVertex3fv(v3);
  glEnd();
}

void subdivide(float *v1, float *v2, float *v3, long depth)
{
  GLfloat v12[3], v23[3], v31[3];
  GLint i;
  if (depth == 0) {
    drawtriangle(v1, v2, v3);
    return;
  }
  for (i = 0; i < 3; i++) {
    v12[i] = v1[i]+v2[i];
    v23[i] = v2[i]+v3[i];
    v31[i] = v3[i]+v1[i];
  }
  normalize(v12);
  normalize(v23);
  normalize(v31);
  subdivide(v1, v12, v31, depth-1);
  subdivide(v2, v23, v12, depth-1);
  subdivide(v3, v31, v23, depth-1);
  subdivide(v12, v23, v31, depth-1);
}

void icosphere(long depth)
{
  #define X .525731112119133606
  #define Z .850650808352039932
  GLfloat vdata[12][3] = {
    {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
    {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
    {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
  };
  GLuint tindices[20][3] = {
//    {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
//    {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
//    {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
//    {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
    {1,4,0}, {4,9,0}, {4,5,9}, {8,5,4}, {1,8,4},
    {1,10,8}, {10,3,8}, {8,3,5}, {3,2,5}, {3,7,2},
    {3,10,7}, {10,6,7}, {6,11,7}, {6,0,11}, {6,1,0},
    {10,1,6}, {11,0,9}, {2,11,9}, {5,2,9}, {11,2,7}
  };

  glEnable(GL_LIGHTING);
  for (int i = 0; i < 20; i++) {
    subdivide(&vdata[tindices[i][0]][0],
    &vdata[tindices[i][1]][0],
    &vdata[tindices[i][2]][0], depth);
  }
  glDisable(GL_LIGHTING);
}

int forward_face(double tilt, double rotate, bool *xp, bool *yp, bool *zp)
{
  *xp = (rotate < 0.0) ? true : false;                         // x positive?
  *yp = ((rotate > 90.0) || (rotate <= -90.0)) ? true : false; // y positive?
  *zp = (tilt > 0.0) ? true : false;                           // z positive?
  if((tilt > 45.0) || (tilt < -45.0)) return 3;               // z facing
  if(((rotate > -45.0) && (rotate <= 45.0))
    || ((rotate <= -135.0) || (rotate > 135.0))) return 2;    // y facing
  return 1;                                                   // x facing
}

void randinit()
{
  srand(time(NULL));
}

double randpm()
{
  return (-0.5 + (1.0 * rand() / RAND_MAX));
}

void HSVtoRGB(HSV *phsv, RGB *prgb)
{
  double hp = phsv->h / 60.0;
  double c = phsv->v * phsv->s;
  double x = c * (1.0 - fabs( fmod(hp, 2.0) -1.0));
  double m = phsv->v - c;

  double rp, gp, bp;
  rp = gp = bp = 0.0;
       if (hp < 1.0) {rp = c  , gp = x  ; bp = 0.0;}
  else if (hp < 2.0) {rp = x  , gp = c  ; bp = 0.0;}
  else if (hp < 3.0) {rp = 0.0, gp = c  ; bp = x  ;}
  else if (hp < 4.0) {rp = 0.0, gp = x  ; bp = c  ;}
  else if (hp < 5.0) {rp = x  , gp = 0.0; bp = c  ;}
  else if (hp < 6.0) {rp = c  , gp = 0.0; bp = x  ;}

  prgb->r = rp + m;
  prgb->g = gp + m;
  prgb->b = bp + m;
}

void popupMessage(QString title, QString message)
{
  QMessageBox::information(0, title, message);
}

void show_help()
{
  QMessageBox::about(0, "",
    "GL_10 v1.0\n"
    "Keypress Commands\n\n"
    "+   increase point size\n"
    "-   decrease point size\n"
    "A   animated rotation run/pause\n"
    "a   axis display toggle\n"
    "B   display boundary (3D) toggle\n"
    "b   bounding box display toggle\n"
    "c   cycle (2D) surface/flat/line\n"
    "c   cycle (3D) full/half/slice/surface/line\n"
    "d   dither points (2D/3D) toggle\n"
    "F   field cycle (1D) Ez/Ez peak/Hy/EzHy\n"
    "F   field cycle (3D) Ex/Ey/Ez/Exyz\n"
    "f    faster simulation update\n"
    "H   help information\n"
    "K   thin (3D) points\n"
    "k   restore (3D) points\n"
    "L   license information\n"
    "o   orthogonal projection toggle\n"
    "Q   increase alpha (2D/3D)\n"
    "q   decrease alpha (2D/3D)\n"
    "R   reset simulation\n"
    "r    run/pause simulation\n"
    "S   stereo display toggle\n"
    "s   slower simulation update\n"
    "T   time display draw\n"
    "t    time field calculation\n"
    "Y   yee cell display toggle (3D)\n"
    "X   exit\n"
  );
}
void show_splash(QWidget *p)
{
  QMessageBox::about(
    p, "",
    "GL_10 v1.01\n"
    "FDTD Modeling/Simulation/Visualization\n"
    "Electro-Magnetic Fields\n\n"
    "Copyright (C) 2005-2013 John Rugis\n\n"
    "This program is free software: you can redistribute it and/or modify it"
    " under the terms of the GNU General Public License as published by the"
    " Free Software Foundation, either version 3 of the License, or (at your"
    " option) any later version.\n\n"
    "This program is distributed in the hope that it will be useful, but"
    " WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY"
    " or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License"
    " for more details.\n\n"
    "http://www.gnu.org/licenses\n\n"
    "rugis@msu.edu\n"
    "After closing this message box, press H for help."
  );
  p->setFocus();     // give focus (back) to primary widget
  p->grabKeyboard(); // make keypress functional even if pointer not over widget
}
