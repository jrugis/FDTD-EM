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

#include "leaplistener.h"

void cLeapListener::onInit(const Controller& controller) {
//  std::cout << "Initialized" << std::endl;
}

void cLeapListener::onConnect(const Controller& controller) {
//  std::cout << "Connected" << std::endl;
}

void cLeapListener::onDisconnect(const Controller& controller) {
//  std::cout << "Disconnected" << std::endl;
}

void cLeapListener::onExit(const Controller& controller) {
//  std::cout << "Exited" << std::endl;
}

void cLeapListener::onFrame(const Controller& controller) {
  // Get the most recent frame and report some basic information
  const Frame frame = controller.frame();
//  std::cout << "Frame id: " << frame.id()
//            << ", timestamp: " << frame.timestamp()
//            << ", hands: " << frame.hands().count()
//            << ", fingers: " << frame.fingers().count()
//            << ", tools: " << frame.tools().count() << std::endl;

  if (!frame.hands().empty()) {
//    // Get the first hand
    const Hand hand = frame.hands()[0];

    // Check if the hand has any fingers
    const FingerList fingers = hand.fingers();
    if (!fingers.empty()) {
      // Calculate the hand's average finger tip position
      Vector avgPos;
      for (int i = 0; i < fingers.count(); ++i) {
        avgPos += fingers[i].tipPosition();
      }
      avgPos /= (float)fingers.count();
//      std::cout << "Hand has " << fingers.count()
//                << " fingers, average finger tip position" << avgPos << std::endl;
    }

    // Get the hand's sphere radius and palm position
//    std::cout << "Hand sphere radius: " << hand.sphereRadius()
//              << " mm, palm position: " << hand.palmPosition() << std::endl;

    // Get the hand's normal vector and direction
    const Vector normal = hand.palmNormal();
    const Vector direction = hand.direction();

    // Calculate the hand's pitch, roll, and yaw angles
//    std::cout << "Hand pitch: " << direction.pitch() * RAD_TO_DEG << " degrees, "
//              << "roll: " << normal.roll() * RAD_TO_DEG << " degrees, "
//              << "yaw: " << direction.yaw() * RAD_TO_DEG << " degrees" << std::endl << std::endl;
  }
}
