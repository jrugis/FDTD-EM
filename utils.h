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

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <QGLWidget>
#include <QString>

typedef struct
{
  double h; // hue [0.0,360.0)
  double s; // saturation [0.0,1.0]
  double v; // lightness [0.0,1.0]
} HSV;

typedef struct
{
  double r; // red [0.0,1.0]
  double g; // green [0.0,1.0]
  double b; // blue [0.0,1.0]
} RGB;

void HSVtoRGB(HSV *phsl, RGB *prgb);

void randinit();
double randpm();
void fatalError(QString message);
void popupMessage(QString title, QString message);
void show_splash(QWidget *p);
void show_help();
int forward_face(double tilt, double rotate, bool *xp, bool *yp, bool *zp);
void icosphere(long depth);

#define BBOX_GEN(bb, lw, r, g, b, x_min, y_min, z_min, x_max, y_max, z_max)\
  bb = glGenLists(1);\
  glNewList(bb, GL_COMPILE);\
    glLineWidth(lw);\
    glBegin(GL_LINES);\
      glColor3d(r, g, b);\
      glVertex3d(x_min, y_min, z_min);\
      glVertex3d(x_max, y_min, z_min);\
      glVertex3d(x_max, y_min, z_min);\
      glVertex3d(x_max, y_max, z_min);\
      glVertex3d(x_max, y_max, z_min);\
      glVertex3d(x_min, y_max, z_min);\
      glVertex3d(x_min, y_max, z_min);\
      glVertex3d(x_min, y_min, z_min);\
      glVertex3d(x_min, y_min, z_max);\
      glVertex3d(x_max, y_min, z_max);\
      glVertex3d(x_max, y_min, z_max);\
      glVertex3d(x_max, y_max, z_max);\
      glVertex3d(x_max, y_max, z_max);\
      glVertex3d(x_min, y_max, z_max);\
      glVertex3d(x_min, y_max, z_max);\
      glVertex3d(x_min, y_min, z_max);\
      glVertex3d(x_min, y_min, z_min);\
      glVertex3d(x_min, y_min, z_max);\
      glVertex3d(x_max, y_min, z_min);\
      glVertex3d(x_max, y_min, z_max);\
      glVertex3d(x_min, y_max, z_min);\
      glVertex3d(x_min, y_max, z_max);\
      glVertex3d(x_max, y_max, z_min);\
      glVertex3d(x_max, y_max, z_max);\
    glEnd();\
  glEndList();

#endif // UTILS_H_INCLUDED
