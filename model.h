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

#ifndef MODEL_H
#define MODEL_H

#include <QGLWidget>

class GLWidget;

class CModel
{
public:
  bool display_bbox;
  bool display_yee;
  bool display_boundary;
  bool dither;
  int cut_type, field_type;
  int skip; // display thinning: skip cells
  double alpha_fac;
  bool xp, yp, zp; // x, y, z positive facing?
  int face;        // the near facing axis: 1=x, 2=y, 3=z

  virtual void reset() {}
  virtual void step() {}
  virtual void draw() const {}
  virtual void get_status(QString &s) const {s = "";}
  virtual void inc_field_type() {}
  virtual void inc_cut_type() {}

protected:
  CModel(GLWidget *parent = NULL);

  GLWidget *p; // parent
  GLuint bbox; // OpenGL disply list

  virtual void set_material() {}
};

#endif // MODEL_H
