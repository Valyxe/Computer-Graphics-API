/********************************************************************************
 * 																				*
 *		Project3.c																*
 *																				*
 *		Created on: Mar 15, 2012												*
 *		Modified on: Apr 8, 2012												*
 *      Author: Allan Simmons													*
 *																				*
 *		This project is a demonstration of myAPI.h. Included are multiple		*
 *		functions that create images that are transformed by the functions		*
 *		built in myAPI.h, including handling of projection and clipping.		*
 *																				*
 ********************************************************************************/


#include <GL/glut.h>
#include <stdio.h>
#include "myAPI.h"

static int WINDOW_WIDTH = 400;
static int WINDOW_HEIGHT = 400;

/********************************************************************************
 * 									init()										*
 *******************************************************************************/
//This function initializes aspects of OpenGL for drawing.
void init(void)
{

	initMat(viewMatrix);
	setView(0, 0, 400, 400);
	initMat(vrcMatrix);
	initMat(projMatrix);
	glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluOrtho2D(-1, 1, -1, 1);
    setVRC(0, 0, 200, 0, 0, 0, 0, 1, 0);
    setProjection(PARALLEL, 0, 0, 200, 100, 100, 150, 250);

/*
	glViewport(0,0,WINDOW_WIDTH, WINDOW_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-50, 50, -50, 50, 150, 250);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 200, 0, 0, 0, 0, 1, 0);
*/
}

/********************************************************************************
 * 									myCube()									*
 *******************************************************************************/
//This function draws a cube using distinct lines.
void myCube()
{
	drawLine(0, 0, 0, 1, 0, 0); //Back face
	drawLine(1, 0, 0, 1, 1, 0);
	drawLine(1, 1, 0, 0, 1, 0);
	drawLine(0, 1, 0, 0, 0, 0);
	drawLine(0, 0, 1, 1, 0, 1); //Front face
	drawLine(1, 0, 1, 1, 1, 1);
	drawLine(1, 1, 1, 0, 1, 1);
	drawLine(0, 1, 1, 0, 0, 1);
	drawLine(0, 0, 0, 0, 0, 1); //Connect front to back
	drawLine(1, 0, 0, 1, 0, 1);
	drawLine(1, 1, 0, 1, 1, 1);
	drawLine(0, 1, 0, 0, 1, 1);
}

/********************************************************************************
 * 									display()									*
 *******************************************************************************/
//This function is the display function that calls draw functions.
void display()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glColor3f(0.0, 0.0, 0.0);

	//Draw cube away from origin
	initMat(curMatrix);
	myScale(20, 40, 20);
	myTranslate(50, -30, -40);
	myRotate(45, 0, 1, 0);
	myCube();

	glFlush();
}

/********************************************************************************
 * 									main()										*
 *******************************************************************************/
//This is the main function of the program. Initializes OpenGL and calls the main loop.
int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGBA|GLUT_DEPTH);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(0, 0);
	glutCreateWindow(argv[0]);

	init();

	glutDisplayFunc(display);

	glutMainLoop();
	return 0;
}
