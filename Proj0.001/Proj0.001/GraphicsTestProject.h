/*
#ifndef GRAPHICSTESTPROJECT_H
#define GRAPHICSTESTPROJECT_H

#include <vector>
#include "atom.h"
#include "vec.h"


// angle of rotation for the camera direction
//float angle = 0.0f;
// actual vector representing the camera's direction
//float lx=0.0f,lz=-1.0f;
// XZ position of the camera
//float x=0.0f, z=5.0f;
// the key states. These variables will be zero
//when no key is being presses
//float deltaAngle = 0.0f;
//float deltaMove = 0;

void changeSize(int w, int h);
void drawSnowMan();
void drawAtom();

void computePos(float deltaMove);

void computeHei(float deltaMoveY);

void computeDir(float deltaAngle);

void renderScene(void);

void pressKey(int key, int xx, int yy);

void releaseKey(int key, int x, int y); 

void processNormalKeys(unsigned char key, int x, int y);

void plotter(int argc, char** argv,std::vector<Atom*>,int);



#endif
*/