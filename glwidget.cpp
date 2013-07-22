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

#include <cstdio>
#include <math.h>

#include "defs.h"
#include "utils.h"
#include "window.h"
#include "model.h"
#include "model1d.h"
#include "model2d.h"
#include "model3d.h"
#include "spaceEH1d.h"
#include "spaceEH2d.h"
#include "spaceEH3d.h"
#include "stopwatch.h"

#include "glwidget.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

// ***********************************************************************
// visual display constants
// ***********************************************************************
#define MSP 50           // DEFAULT MODEL STEP PERIOD (ms)
#define SPCF 1.2         // STEP PERIOD CHANGE FACTOR
#define ARR 20.0         // ANIMATE_REFRESH_RATE
#define AP (1000.0/ARR)  // ANIMATE_PERIOD (ms)
#define RAA (0.005)      // ROTATE_ANIMATION_ANGLE per ms
#define RDS 0.2          // ROTATE_DRAG_SCALE

#define INIT_ROTATE (0.0)
#define ROTATE_MIN (-180.0)
#define ROTATE_MAX (180.0)

#define INIT_TILT (0.0)
#define TILT_MIN (-90.0)
#define TILT_MAX (90.0)

#define D13_1440x900
//#define D27_1920x1280

// **** 27 inch display (1920x1280)
#ifdef D27_1920x1280
  #define ED 0.8    // EYE_DISTANCE (meters)
  #define SH 0.285   // SCREEN_HEIGHT
  #define OH 0.175   // OBJECT_HEIGHT
#endif

// **** 13 inch display (1440x900 MacBook Air)
#ifdef D13_1440x900
  #define ED 0.6    // EYE_DISTANCE (meters)
  #define SH 0.176  // SCREEN_HEIGHT
  #define OH 0.1    // OBJECT_HEIGHT
#endif

#define TT   ((SH/2.0)/ED)  // TAN_THETA
#define FRN  (ED*0.1)       // FRUSTUM_NEAR
#define FRF  (ED*1.9)       // FRUSTUM_FAR
#define FRL  (-1.0*TT*FRN)  // FRUSTUM_LEFT
#define FRR  (TT*FRN)       // FRUSTUM_RIGHT
#define FRB  (-1.0*TT*FRN)  // FRUSTUM_BOTTOM
#define FRT  (TT*FRN)       // FRUSTUM_TOP
#define ON   FRN            // ORTHO_NEAR
#define OF   FRF            // ORTHO_FAR
#define OL   (-1.0*SH/2.0)  // ORTHO_LEFT
#define OR   (SH/2.0)       // ORTHO_RIGHT
#define OB   (-1.0*SH/2.0)  // ORTHO_BOTTOM
#define OT   (SH/2.0)       // ORTHO_TOP

#define SCALE_BUMP 0.0005
#define SCALE_MIN (0.1*OH)
#define SCALE_MAX (5.0*OH)

// stereo parameters
//#define ES 0.07   // EYE_SEPARATION
#define SSFP (0.07)    // STEREO_SEP_FACTOR_P
#define SSFN (-0.07)   // STEREO_SEP_FACTOR_N
#define SSF  (0.55)     // STEREO_SCALE_FACTOR

#ifdef USE_LEAP
#include "Leap.h"

#define LEAP_POINTER_MIN 30.0
#define LEAP_RDS (5.0 * RDS) // LEAP ROTATE_DRAG_SCALE
#define LEAP_SCALE_BUMP (15.0 * SCALE_BUMP)
#define LEAP_POINT_SCALE 0.1
#define LEAP_ALPHA_SCALE 0.05
#define SWIPE_MAG_MIN 100.0     // mm
#define SWIPE_SPEED_MIN 2000.0
#define CIRCLE_RADIUS_MIN 75.0 // mm
#endif

GLWidget::GLWidget(QWidget *parent)
  : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
  setFocusPolicy(Qt::ClickFocus); // accept key presses

  model = NULL;
  model_dim = 3;
  animate = false;
  run = false;
#ifdef USE_LEAP
  leap_drag = leap_drag_start = false;
  leap_mode = 3;
#endif
  time_all = time_field = time_paint = false;
  display_ortho = false;
  display_axis = true;
  display_stereo = false;

  scale_factor = OH;
  point_size = 1.0;
  rotate_angle = INIT_ROTATE;
  tilt_angle = INIT_TILT;
  step_period = MSP;

  animation_timer = new QTimer(this);
  connect(animation_timer, SIGNAL(timeout()), this, SLOT(animation_tick()));
  calcs_timer = new QTimer(this);
  connect(calcs_timer, SIGNAL(timeout()), this, SLOT(calcs_tick()));
  EH_timer = new CStopwatch("calc:");
  Paint_timer = new CStopwatch("draw:");
}

GLWidget::~GLWidget()
{
}

#ifdef USE_LEAP
void GLWidget::onInit(const Leap::Controller& controller) {
  //controller.enableGesture(Leap::Gesture::TYPE_KEY_TAP);
  controller.enableGesture(Leap::Gesture::TYPE_SWIPE);
  controller.enableGesture(Leap::Gesture::TYPE_CIRCLE);
}
void GLWidget::onConnect(const Leap::Controller& controller) {
}
void GLWidget::onDisconnect(const Leap::Controller& controller) {
}
void GLWidget::onExit(const Leap::Controller& controller) {
}

void GLWidget::onFrame(const Leap::Controller& controller) {
  const Leap::Frame frame = controller.frame(); // most recent frame

  if(!leap_drag) {
    for(int i = 0; i < frame.gestures().count(); i++) {
      const Leap::Gesture gesture = frame.gestures()[i];
      if(gesture.type() == Leap::Gesture::TYPE_SWIPE) {
        Leap::SwipeGesture swipe(gesture);
        if(fabs(swipe.direction().z) > 0.2) continue;
        //if(swipe.speed() < SWIPE_SPEED_MIN) continue;
        //if((swipe.startPosition() - swipe.position()).magnitude() < SWIPE_MAG_MIN) continue;
        float angle = swipe.direction().roll();
        leap_mode = 1.0 + (6.0 * (1.0 + angle / -PI));
        /*
        if(angle <= -0.917 * PI) leap_mode = 12; // 12 o'clock
        else if(angle <= -0.75 * PI) leap_mode = 11; // 11
        else if(angle <= -0.583 * PI) leap_mode = 10; // 10
        else if(angle <= -0.417 * PI) leap_mode = 9; // 9
        else if(angle <= -0.25 * PI) leap_mode = 8;   // 8
        else if(angle <= -0.0833 * PI) leap_mode = 7; // 7
        else if(angle <= 0.0833 * PI) leap_mode = 6;
        else if(angle <= 0.25 * PI) leap_mode = 5;
        else if(angle <= 0.417 * PI) leap_mode = 4;
        else if(angle <= 0.583 * PI) leap_mode = 3;
        else if(angle <= 0.75 * PI) leap_mode = 2;
        else if(angle <= 0.917 * PI) leap_mode = 1;
        else leap_mode = 12; // 12
        */
        update_status();
      }
      else if(gesture.type() == Leap::Gesture::TYPE_CIRCLE) {
        Leap::CircleGesture circle(gesture);
        if(circle.radius() < CIRCLE_RADIUS_MIN) continue;
        if((circle.normal().z < 0.0) & !run) {  // clockwise
          run = true;
          calcs_timer->start(step_period);
          update_status();
        }
        if((circle.normal().z > 0.0) & run) {  // counter-clockwise
          run = false;
          update_status();
        }
      }
    }
    return;
  }

  if(frame.pointables().count() != 1) return;
  const Leap::Pointable pointer = frame.pointables()[0];
  if(pointer.length() < LEAP_POINTER_MIN) return;
  Leap::Vector fpos = pointer.tipPosition();
  if(leap_drag_start) {
    flastPos = fpos;
    leap_drag_start = false;
    return;
  }
  Leap::Vector fmotion = fpos - flastPos;
  flastPos = fpos;

  switch(leap_mode) {
  case 3: // spin
    tilt_angle -= fmotion.y * LEAP_RDS;
    if(tilt_angle > TILT_MAX) tilt_angle = TILT_MAX;
    if(tilt_angle < TILT_MIN) tilt_angle = TILT_MIN;

    rotate_angle += fmotion.x * LEAP_RDS;
    if(rotate_angle > ROTATE_MAX) rotate_angle -= 360.0;
    if(rotate_angle < ROTATE_MIN) rotate_angle += 360.0;

    scale_factor *= 1.0 + (fmotion.z * LEAP_SCALE_BUMP);
    if(scale_factor < SCALE_MIN) scale_factor = SCALE_MIN;
    if(scale_factor > SCALE_MAX) scale_factor = SCALE_MAX;
    break;

  case 12: // points
    point_size += LEAP_POINT_SCALE * fmotion.y;
    if(point_size < 1.0) point_size = 1.0;
    model->alpha_fac *= 1.0 + LEAP_ALPHA_SCALE * fmotion.x;
    if(!run && !animate) updateGL();
    break;
  }

  if(!run && !animate) updateGL();
}

void GLWidget::get_leap_status(QString &s)
{
  s.append(", LEAP");
  if(leap_mode == 3) s.append("[spin]");
  else if(leap_mode == 12) s.append("[points]");
  else s.append("[" + QString::number(leap_mode) + "]");
}
#endif

void GLWidget::keyPressEvent(QKeyEvent *e)
{
  if(e->text() == "1") {
    model_dim = 1;
    reset_model();
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "2") {
    model_dim = 2;
    reset_model();
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "3") {
    model_dim = 3;
    reset_model();
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "A") {
    animate ^= 1;
    if(animate) {
      animation_timer->start(AP);
      elapsed_time.restart();
    }
    else animation_timer->stop();
  }
  else if(e->text() == "a") {
    display_axis ^= 1;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "B") {
    model->display_boundary ^= 1;
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "b") {
    model->display_bbox ^= 1;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "c") {
    model->inc_cut_type();
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "d") {
    model->dither ^= 1;
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "F") {
    model->inc_field_type();
    update_status();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "f") {
    step_period /= SPCF;
    if(run) calcs_timer->start(step_period);
  }
  else if(e->text() == "H") {
    show_help();
  }
  else if(e->text() == "K") {
    model->skip++;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "k") {
    if(model->skip > 1) {
      model->skip--;
      if(!run && !animate) updateGL();
    }
  }
  else if(e->text() == "L") {
    show_splash(this);
  }
  else if(e->text() == "Q") {
    model->alpha_fac *= 1.1;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "q") {
    model->alpha_fac /= 1.1;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "o") {
    display_ortho ^= 1;
    update_status();
    reset_projection();
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "r") toggle_run();
  else if(e->text() == "R") {
    model->reset();
    point_size = 1.0;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "s") {
    step_period *= SPCF;
    if(run) calcs_timer->start(step_period);
  }
  else if(e->text() == "S") {
    display_stereo ^= 1;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "t") {
    time_all ^= 1;
    update_status();
  }
  else if(e->text() == "Y") {
    model->display_yee ^= 1;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "+") {
    point_size += 0.5;
    if(!run && !animate) updateGL();
  }
  else if(e->text() == "-") {
    if(point_size > 1.0) {
      point_size -= 0.5;
      if(!run && !animate) updateGL();
    }
  }
  else if(e->text() == "X") {
    parentWidget()->close();
  }
#ifdef USE_LEAP
  else if(e->key() == Qt::Key_Control) {
    leap_drag = leap_drag_start = true;
  }
#endif
  else QGLWidget::keyPressEvent(e);
}

void GLWidget::keyReleaseEvent(QKeyEvent *e)
{
#ifdef USE_LEAP
  if(e->key() == Qt::Key_Control) {
    leap_drag = leap_drag_start = false;
  }
  else
#endif
  QGLWidget::keyReleaseEvent(e);
}

void GLWidget::toggle_run()
{
  run ^= 1;
  update_status();
  if(run) calcs_timer->start(step_period);
}

void GLWidget::initializeGL()
{
  glEnable(GL_DOUBLEBUFFER);
  glEnable(GL_CULL_FACE);
  glShadeModel(GL_SMOOTH);
  //glShadeModel(GL_FLAT);

  GLfloat light_ambient[] = {0.5, 0.5, 0.5, 1.0};
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambient);

  GLfloat light_diffuse[] = {0.9, 0.9, 0.9, 1.0};
  GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_position[] = {0.0, 0.0, 5.0, 1.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.3);
  glEnable(GL_LIGHT0);

  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  //glEnable(GL_MULTISAMPLE);
  //glEnable(GL_MULTISAMPLE_ARB);
  //glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable (GL_SAMPLE_ALPHA_TO_COVERAGE); // includes depth buffer in alpha blending !!!
  //glAlphaFunc(GL_GREATER, 0.5);
  //glEnable(GL_ALPHA_TEST);

  glPointSize(point_size);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH, GL_NICEST);
  //glHint(GL_LINE_SMOOTH, GL_FASTEST);
  //glHint(GL_LINE_SMOOTH, GL_DONT_CARE);

  iwidth = cwidth = width(); // initial pixel dimension of window
  iheight = cheight = height();
  reset_projection();
  glMatrixMode(GL_MODELVIEW); // setup MODEL matrix
  glLoadIdentity();

  axes = glGenLists(1);
  glNewList(axes, GL_COMPILE);
    glLineWidth(3.0);
    glBegin(GL_LINES);
      glColor3d(MIDRED);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.5, 0.0, 0.0);
      glColor3d(MIDGREEN);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.0, 0.5, 0.0);
      glColor3d(MIDBLUE);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.0, 0.0, 0.5);
    glEnd();
  glEndList();

  reset_model(); // here because relies on OpenGL
  update_status();
  show_splash(this);
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
  int dx = event->x() - lastPos.x();
  int dy = event->y() - lastPos.y();

  if(event->buttons() & Qt::LeftButton) {
    rotate_angle += dx * RDS;
    if(rotate_angle > ROTATE_MAX) rotate_angle -= 360.0;
    if(rotate_angle < ROTATE_MIN) rotate_angle += 360.0;

    tilt_angle += dy * RDS;
    if(tilt_angle > TILT_MAX) tilt_angle = TILT_MAX;
    if(tilt_angle < TILT_MIN) tilt_angle = TILT_MIN;
    updateGL();
  }
  lastPos = event->pos();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
  lastPos = event->pos();
  if(event->buttons() & Qt::LeftButton) {
    if(event->pos().x() < 20) toggle_run();
  }
}

void GLWidget::draw()
{
  glPointSize(point_size);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glRotated(tilt_angle, 1,0,0);
  glRotated(rotate_angle, 0,0,1);
  glScaled(scale_factor, scale_factor, scale_factor);

  model->face = forward_face(tilt_angle, rotate_angle, &(model->xp), &(model->yp), &(model->zp));
  if(display_axis) glCallList(axes);
  model->draw();
  glPopMatrix();

 //*******************************************
 // little red box
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho (0, this->size().width(), this->size().height(), 0, 0, 1);
  //glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glColor4f(1.0, 0.0, 0.0, 0.5);
  glBegin(GL_QUADS);
  glVertex2f(0, 0); glVertex2f(0, 10); glVertex2f(10, 10); glVertex2f(10, 0);
  glEnd();

  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  //*******************************************
}

void GLWidget::paintGL()
{
  if(time_all) time_paint = true;
  if(time_paint) Paint_timer->start();  // time the paint?
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(display_stereo) {
    glMatrixMode(GL_PROJECTION);
    // right eye, left side
    glPushMatrix();
    glTranslatef(SSFN, 0, 0);
    glScalef(SSF, SSF, SSF);
    draw();
    glPopMatrix();
    // left eye, right side
    glPushMatrix();
    glTranslatef(SSFP, 0, 0);
    glScalef(SSF, SSF, SSF);
    draw();
    glPopMatrix();
  } else draw();

  if(time_paint) { // display the elapsed time
    Paint_timer->stop();
    time_paint = false;
    update_status();
  }
}

void GLWidget::resizeGL(int width, int height)
{
  int side = qMax(width, height);
  glViewport((width - side) / 2, (height - side) / 2, side, side);
  cwidth = width;
  cheight = height;
  reset_projection();
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
  scale_factor *= 1.0 + (SCALE_BUMP * event->delta());
  if(scale_factor < SCALE_MIN) scale_factor = SCALE_MIN;
  if(scale_factor > SCALE_MAX) scale_factor = SCALE_MAX;
  updateGL();
}

void GLWidget::animation_tick()
{
  rotate_angle += RAA * elapsed_time.elapsed();
  elapsed_time.restart();
  if(rotate_angle > ROTATE_MAX) rotate_angle -= 360.0;
  updateGL();
}

void GLWidget::calcs_tick()
{
  if(time_all) time_field = true;
  if(time_field) EH_timer->start();  // time the EH field update?
  model->step();
  if(time_field) { // display the elapsed time
    EH_timer->stop();
    time_field = false;
    update_status();
  }
  if(!run) calcs_timer->stop();
  if(!animate) updateGL();
}

void GLWidget::reset_projection()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if(display_ortho) glOrtho(OL, OR, OB, OT, ON, OF);
  else glFrustum(FRL, FRR, FRB, FRT, FRN, FRF);
  glTranslated(0.0, 0.0, -1.0 * ED);
  glRotated(-90.0, 1,0,0);
  glScalef(1.0 / WSF, 1.0 / WSF, 1.0 / WSF);
  //float s = qMax(iwidth / cwidth / WSF, iheight / cheight / WSF);
  //glScalef(s, s, 1.0 / WSF);
  //float sw = cwidth / iwidth / WSF;
  //float sh = cheight / iheight / WSF;
  //glScalef(sw, sh, 1.0 / WSF);
}

void GLWidget::reset_model()
{
  if(model != NULL) delete model;
  switch(model_dim) {
  case 1:
    glClearColor(BLACK, 1.0);
    model = new CModel1D(this);
    break;
  case 2:
    glClearColor(BLACK, 1.0);
    model = new CModel2D(this);
    break;
  case 3:
    glClearColor(LIGHTGREY, 1.0);
    model = new CModel3D(this);
    break;
  }
}

void GLWidget::update_status() {  // display status in title bar
  QString s = TITLE;
  s.append(": ");
  model->get_status(s);
  if(display_ortho) s.append(", ortho");
  if(run) s.append(", running");
  if(time_all) {
    s.append(", " + EH_timer->message);
    s.append(", " + Paint_timer->message);
  }
#ifdef USE_LEAP
  get_leap_status(s);
#endif
  //main_window->setWindowTitle(s);
  parentWidget()->setWindowTitle(s);
}

