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

#ifndef MODEL2D_H
#define MODEL2D_H

#include "model.h"

class CSpaceEH2d;
class CAbc2o2d;
class CTfsf2d;

class CModel2D : public CModel
{
public:
  CModel2D(GLWidget *parent = 0);

  void reset();
  void step();
  void draw() const;
  void get_status(QString &s) const;
  void inc_field_type();
  void inc_cut_type();

private:
  CSpaceEH2d *space2d; // 2d space
  CAbc2o2d *abc2o2d; // second order abc
  CTfsf2d *tfsf2d; // tfsf in 2d space

  size_t time_step;  // time step
  void set_material();
};

#endif // MODEL2D_H
