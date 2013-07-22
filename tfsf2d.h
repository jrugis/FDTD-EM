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

#ifndef TFSF2D_H
#define TFSF2D_H

class Ccell1d;
class Ccell2d;
class CSpaceEH1d;
class CSpaceEH2d;
class CAbc2o1d;

class CTfsf2d
{
public:
  CTfsf2d(CSpaceEH2d *s, size_t sB, size_t sD);
  void reset();
  void updateA();     // before source signal update
  void updateB();     // after source signal update
  double *inp, *inpm1; // source signal inputs
private:
  CSpaceEH1d *a;     // line wave source: 1D auxillary space
  Ccell2d *c;        // the 2D model space
  CAbc2o1d *abc2o1d; // second order abc
  size_t sx, sy;     // size of 2D model space
  size_t sb;         // size of tfsf boundary in 2D model (per edge)
  size_t sd;         // size 1D aux space decay regions
  size_t sa;         // total size of 1D aux space
};

#endif // TFSF2D_H
