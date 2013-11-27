#include <stdlib.h>
#include <math.h>
#include <vector>
#include "atom.h"
#include "vec.h"


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// angle of rotation for the camera direction
float angle = 0.0f;
// actual vector representing the camera's direction
float lx=0.0f,lz=-1.0f;
// XZ position of the camera
float x=0.0f, z=30.0f;
// Y position of the camera
float y=1.0f;
// the key states. These variables will be zero
//when no key is being presses
float deltaAngle = 0.0f;
float deltaMove = 0;
float deltaMoveY = 0;
float deltaMoveZ = 0;

std::vector<Atom*> list_of_atoms;
int number_of_atoms;

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;
	float ratio =  w * 1.0f / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 500.0f);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

void drawAtom() {
	glColor3f(0.0f, 1.0f, 0.0f); // Make the atoms green
	glutSolidSphere(0.75f,20,20); // Make the atoms spheres
	
}


void computePos(float deltaMove) {

	x += deltaMove * lx * 0.1f;
	z += deltaMove * lz * 0.1f;
}

void computeHei(float deltaMoveY) {

	y += deltaMoveY * 0.1f;
}

void computeSid(float deltaMoveZ) {
	z += deltaMoveZ*0.1f;
}

void computeDir(float deltaAngle) {

	angle += deltaAngle;
	lx = sin(angle);
	lz = -cos(angle);
}

void renderScene(void) {

	if (deltaMove)
		computePos(deltaMove);
	if(deltaMoveY)
		computeHei(deltaMoveY);
	if(deltaMoveZ)
		computeSid(deltaMoveZ);
	if (deltaAngle)
		computeDir(deltaAngle);

	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();
	// Set the camera
	gluLookAt(	x, y, z,
				x+lx, y,  z+lz,
				0.0f, 1.0f,  0.0f);

	// Draw atoms
	for(unsigned int k=0; k < list_of_atoms.size()/number_of_atoms;k++){
		for(int i = 0; i < number_of_atoms;i++){
		glPushMatrix();
		glTranslatef(list_of_atoms[i]->get_position().getX(),list_of_atoms[i]->get_position().getY(),list_of_atoms[i]->get_position().getZ());
		drawAtom();
		glPopMatrix();
		}
		if(k>0){
			system("pause");
		}
	}
	glutSwapBuffers();
}


void pressKey(int key, int xx, int yy) { // HUR RÖRA SIG I Y???

	switch (key) {
		case GLUT_KEY_LEFT : deltaAngle = -0.01f; break;
		case GLUT_KEY_RIGHT : deltaAngle = 0.01f; break;
		case GLUT_KEY_UP : deltaMove = 0.5f; break;
		case GLUT_KEY_DOWN : deltaMove = -0.5f; break;
		case GLUT_KEY_PAGE_UP : deltaMoveY = 0.5f; break;
		case GLUT_KEY_PAGE_DOWN : deltaMoveY = -0.5f; break;
		case GLUT_KEY_INSERT : deltaMoveZ = -0.5f; break;
		case GLUT_KEY_HOME : deltaMoveZ = 0.5f; break;
	}
}


void releaseKey(int key, int xxx, int yyy) {
	switch (key) {
		case GLUT_KEY_LEFT :
		case GLUT_KEY_RIGHT : deltaAngle = 0.0f;break;
		case GLUT_KEY_UP :
		case GLUT_KEY_DOWN : deltaMove = 0;break;
		case GLUT_KEY_PAGE_UP :
		case GLUT_KEY_PAGE_DOWN : deltaMoveY = 0; break;
		case GLUT_KEY_INSERT :
		case GLUT_KEY_HOME : deltaMoveZ = 0; break;
	}
}

void plotter(int argc, char** argv,std::vector<Atom*> incoming_list_of_atoms,int incoming_number_of_atoms) {

	list_of_atoms=incoming_list_of_atoms;
	number_of_atoms=incoming_number_of_atoms;
	// init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(320,320);
	glutCreateWindow("Atoms");
	glClearColor(1,1,1,1); // Make the background white
	// register callbacks
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);
	glutSpecialFunc(pressKey);
	// here are the new entries
	glutIgnoreKeyRepeat(1);
	glutSpecialUpFunc(releaseKey);
	// OpenGL init
	glEnable(GL_DEPTH_TEST);
	// enter GLUT event processing cycle
	glutMainLoop();
}