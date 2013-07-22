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

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QtGui>
#include <QWidget>
#include <QtOpenGL>
#include <QGLWidget>
#include <QDialog>
#include <QString>
#include <QTime>

#include "defs.h"
#ifdef USE_LEAP
#include "Leap.h"
#endif

class CModel;
class Window;
class CStopwatch;

#ifdef USE_LEAP
class GLWidget : public QGLWidget, public Leap::Listener
#else
class GLWidget : public QGLWidget
#endif
{
  Q_OBJECT

public:
  GLWidget(QWidget *parent);
  ~GLWidget();
  bool run;

public slots:
  void animation_tick();
  void calcs_tick();
signals:
  void timeout();

protected:
  void initializeGL();
  void keyPressEvent(QKeyEvent *e);
  void keyReleaseEvent(QKeyEvent *e);
  void mouseMoveEvent(QMouseEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void paintGL();
  void resizeGL(int width, int height);
  void wheelEvent(QWheelEvent *event);

#ifdef USE_LEAP
  virtual void onInit(const Leap::Controller&);
  virtual void onConnect(const Leap::Controller&);
  virtual void onDisconnect(const Leap::Controller&);
  virtual void onExit(const Leap::Controller&);
  virtual void onFrame(const Leap::Controller&);
#endif

private:
#ifdef USE_LEAP
  bool leap_drag, leap_drag_start;
  int leap_mode;
  Leap::Vector flastPos;
  void get_leap_status(QString &s);
#endif

  QTimer *animation_timer, *calcs_timer;
  QTime elapsed_time;
  QPoint lastPos;
  GLdouble scale_factor;
  GLdouble point_size;
  GLdouble rotate_angle; // +-180 degrees (wraps)
  GLdouble tilt_angle;   // +-90 degrees (locks)
  GLuint axes;
  //GLuint iwidth, iheight; // initial pixel dimension of window
  //GLuint cwidth, cheight; // current pixel dimension of window
  float iwidth, iheight; // initial pixel dimension of window
  float cwidth, cheight; // current pixel dimension of window

  QString title;
  int model_dim;
  bool animate;
  bool display_ortho, display_axis, display_stereo;
  bool time_all, time_field, time_paint;
  double step_period;
  CModel *model;
  CStopwatch *EH_timer, *Paint_timer;

  void draw();
  void reset_model();
  void reset_projection();
  void update_status();

  void toggle_run();

};

#endif
