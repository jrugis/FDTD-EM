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

#include <QApplication>
#include <QDesktopWidget>
#include <QGLFormat>

#include "utils.h"
#include "defs.h"
#include "window.h"

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);
  if (!QGLFormat::hasOpenGL()) fatalError("This system does not support OpenGL.");
  Window window;
  int h = QApplication::desktop()->height() * WSF;
  //int w = QApplication::desktop()->width() * WSF;
  window.resize(h, h);
  window.show();
  //window.showMaximized();
  window.setMinimumSize(window.size());
  window.setMaximumSize(window.size());
  return app.exec();
}


