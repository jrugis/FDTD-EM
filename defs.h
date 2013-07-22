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

#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

/*
   dt = dx / c
   l = c / f
   where: dx (m)
          dt (s)
          f (1/s)
          l (m)
          c = 3x10^8 (m/s)


   example:
     SIZEX = 300
     SPACE = 0.2 (m)
     f = 2.4x10^9

     dx = SPACE / SIZEX = 6.66X10^-5
     dt = 2.2x10^-13
     l = 0.125   (~= SPACE ?)
     l / dx =    (>20, ?)
*/

#define TITLE "GL_10"
#define USE_LEAP

// ***********************************************************************
// visual display constants
// ***********************************************************************
#define WSF 0.93 // WINDOW_SIZE_FACTOR

// ***********************************************************************
// modelling defs
// ***********************************************************************
#define D22    // field update method
//#define D24

// ***********************************************************************
// electrical constants
// ***********************************************************************
#define IMP0 377.0
#define DTDS2D (1.0 / sqrt(2.0))
#define DTDS3D (1.0 / sqrt(3.0))

// ***********************************************************************
// math constants
// ***********************************************************************
#define PI 3.14159265358979323846

// ***********************************************************************
// OpenGL constants
// ***********************************************************************
#define WHITE 1.0, 1.0, 1.0
#define LIGHTGREY 0.8, 0.8, 0.8
#define MIDGREY 0.5, 0.5, 0.5
#define DARKGREY 0.3, 0.3, 0.3
#define BLACK 0.0, 0.0, 0.0
#define BRIGHTRED 1.0, 0.0, 0.0
#define BRIGHTGREEN 0.0, 1.0, 0.0
#define BRIGHTBLUE 0.0, 0.0, 1.0
#define MIDRED 0.6, 0.0, 0.0
#define MIDGREEN 0.0, 0.6, 0.0
#define MIDBLUE 0.0, 0.0, 0.6
#define BRIGHTCYAN 0.0, 1.0, 1.0
#define BRIGHTMAGENTA 1.0, 0.0, 1.0
#define BRIGHTYELLOW 1.0, 1.0, 0.0
#define MIDCYAN 0.0, 0.4, 0.4
#define MIDMAGENTA 0.4, 0.0, 0.4
#define MIDYELLOW 0.4, 0.4, 0.0

// some colors as GLdouble arrays
//const GLdouble black[] = {0.0, 0.0, 0.0, 1.0};
//const GLdouble white[] = {1.0, 1.0, 1.0, 1.0};
//const GLdouble red[] = {0.9, 0.0, 0.0, 1.0};
//const GLdouble green[] = {0.0, 0.9, 0.0, 1.0};
//const GLdouble blue[] = {0.0, 0.0, 0.9, 1.0};
//const GLdouble yellow[] = {0.9, 0.9, 0.0, 1.0};
//const GLdouble pink[] = {1.0, 0.73, 1.0, 1.0};
//const GLdouble orange[] = {0.8, 0.4, 0.0, 1.0};

#endif // DEFS_H_INCLUDED
