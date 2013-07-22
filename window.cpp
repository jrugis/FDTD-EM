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

#include <QtGui>
#include <Qt>

#include "defs.h"
#include "glwidget.h"

#ifdef USE_LEAP
#include "Leap.h"
#endif

#include "window.h"

Window::Window()
{
  setFocusPolicy(Qt::ClickFocus); // accept key-presses
  glWidget = new GLWidget(this);
#ifdef USE_LEAP
  controller = new Leap::Controller;
  controller->addListener(*glWidget);
#endif
  QHBoxLayout *mainLayout = new QHBoxLayout;
  mainLayout->setContentsMargins(QMargins(5, 5, 5, 5));
  mainLayout->addWidget(glWidget);
  setLayout(mainLayout);
  setWindowTitle(tr("GL_10"));
}

Window::~Window()
{
#ifdef USE_LEAP
  controller->removeListener(*glWidget);
#endif
}
