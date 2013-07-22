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

#ifndef TFSF3D_H
#define TFSF3D_H

class Ccell1d;
class Ccell3d;
class CSpaceEH1d;
class CSpaceEH3d;
class CAbc2o1d;

class CTfsf3d
{
public:
  CTfsf3d(CSpaceEH3d *s, size_t sB, size_t sD);
  void reset();
  void updateA();      // before source signal update
  void updateB();      // after source signal update
  double *inp, *inpm1; // source signal inputs
  size_t sb;   // size of tfsf boundary in 3D model (per face)
private:
  CSpaceEH1d *a;       // plane wave source: 1D auxillary space
  Ccell3d *c;          // the 3D model space
  CAbc2o1d *abc2o1d;   // second order abc
  size_t sx, sy, sz;   // size of 3D model space
  size_t sxy;
  size_t sd;   // size 1D aux space decay region
  size_t sa;   // total size of 1D aux space
};

#endif // TFSF3D_H
