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

#ifndef MODEL3D_H
#define MODEL3D_H

#include "model.h"

class CSpaceEH3d;
class CAbc1o3d;
class CTfsf3d;

class CModel3D : public CModel
{
public:
  CModel3D(GLWidget *parent = 0);

  void reset();
  void step();
  void draw() const;
  void get_status(QString &s) const;
  void inc_field_type();
  void inc_cut_type();

private:
  CSpaceEH3d *space3d; // 3d space
  CAbc1o3d *abc1o3d; // first order abc
  CTfsf3d *tfsf3d; // tfsf in 3d space
  GLuint objects;

  size_t time_step;  // time step
  void set_material();
};

#endif // MODEL3D_H
