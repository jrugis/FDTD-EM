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

#include "utils.h"
#include "glwidget.h"

#include "model.h"

CModel::CModel(GLWidget *parent)
{
  p = parent;
  display_bbox = true;
  dither = false;
  display_yee = false;
  display_boundary = false;
  alpha_fac = 1.0;
  skip = 1;
  cut_type = field_type = 0;

  // create data bounding box
  BBOX_GEN(bbox, 3.0, 0.3, 0.3, 0.3, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5)
}
