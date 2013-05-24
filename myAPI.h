/********************************************************************************\
 *																				*
 *		myAPI.h																	*
 *																				*
 *		Created on: Feb 25, 2012												*
 *		Modified on: Apr 1, 2012												*
 *      Author: Allan Simmons													*
 *      																		*
 *      This API includes functions that handle transformation of points for	*
 *      OpenGL. The transforms and points are handled as matrices, and the		*
 *      actual transformations are handled with the multiplication of transform	*
 *      matrices by the point matrices.											*
 *      																		*
 ********************************************************************************/

#include <math.h>
#include <string.h>
#include <GL/glut.h>

#ifndef MYAPI_H_
#define MYAPI_H_

typedef double myVertex[4];
typedef double myMat[4][4];

void initMat(myMat);
int clipLine(myVertex, myVertex, double, double, double, double, myVertex, myVertex);

myMat curMatrix;
myMat vrcMatrix;
myMat projMatrix;
myMat viewMatrix;

double vrp[3];
double vup[3];

int TYPE;
enum TYPE {PARALLEL, PERSPECTIVE};
/********************************************************************************
 * 									multMatrix4x4()								*
 *******************************************************************************/
//This function multiplies two matrices, mat1 and mat2, storing the results in mat1.
void multMatrix4x4(myMat mat1, myMat mat2)
{
	myMat tempMat;
	int x, y;
	for(x = 0; x < 4; x++)
	{
		for(y = 0; y < 4; y++)
		{
			tempMat[x][y] = 0;
		}
	}

	int i, j, k;
	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			for(k = 0; k < 4; k++)
			{
				tempMat[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}

	memcpy(mat1, tempMat, 16*sizeof(double));
	//printf("%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n\n", mat1[0][0], mat1[1][0], mat1[2][0], mat1[3][0], mat1[0][1], mat1[1][1], mat1[2][1], mat1[3][1], mat1[0][2], mat1[1][2], mat1[2][2], mat1[3][2], mat1[0][3], mat1[1][3], mat1[2][3], mat1[3][3]);
}

/********************************************************************************
 * 									transformVertex()							*
 *******************************************************************************/
//This function transforms a point by multiplying its matrix by the transform matrix.
double* transformVertex(myMat matrix, myVertex vertex)
{
	double* newVector = malloc(sizeof(myVertex));

	int x;
	for(x = 0; x < 4; x++)
	{
		newVector[x] = 0;
	}

	int i, j;
	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			newVector[i] += matrix[j][i] * vertex[j];
		}
	}

	return newVector;
}

/********************************************************************************
 * 									initMat()									*
 *******************************************************************************/
//This function initializes a matrix to the identity matrix.
void initMat(myMat matrix)
{
	int i, j;
	for(i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			if(i == j)
			{
				matrix[i][j] = 1;
			}
			else
			{
				matrix[i][j] = 0;
			}
		}
	}
}

/********************************************************************************
 * 									myTranslate()								*
 *******************************************************************************/
//This function creates a translate transform matrix.
void myTranslate(double dx, double dy, double dz)
{
	myMat transMat;
	initMat(transMat);

	transMat[3][0] = dx;
	transMat[3][1] = dy;
	transMat[3][2] = dz;

	multMatrix4x4(curMatrix, transMat);
}

/********************************************************************************
 * 									myScale()									*
 *******************************************************************************/
//This function creates a scale transform matrix.
void myScale(double sx, double sy, double sz)
{
	myMat scaleMat;
	initMat(scaleMat);

	scaleMat[0][0] = sx;
	scaleMat[1][1] = sy;
	scaleMat[2][2] = sz;

	multMatrix4x4(curMatrix, scaleMat);
}

/********************************************************************************
 * 									myRotate()									*
 *******************************************************************************/
//This function creates a rotate transform matrix based on quaternions.
void myRotate(double theta, double x, double y, double z)
{
	myMat rotateMat;
	initMat(rotateMat);

	int vector[] = {x, y, z};
	int vectorMag = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	int i;
	for(i = 0; i < 3; i++)
	{
		vector[i] /= vectorMag;
	}

	rotateMat[0][0] = pow(vector[0], 2)*(1 - cos(theta)) + cos(theta);				rotateMat[1][0] = vector[0]*vector[1]*(1 - cos(theta)) - vector[2]*sin(theta);	rotateMat[2][0] = vector[0]*vector[2]*(1 - cos(theta)) + vector[1]*sin(theta);
	rotateMat[0][1] = vector[1]*vector[0]*(1 - cos(theta)) + vector[2]*sin(theta);	rotateMat[1][1] = pow(vector[1], 2)*(1 - cos(theta)) + cos(theta);				rotateMat[2][1] = vector[1]*vector[2]*(1 - cos(theta)) - vector[0]*sin(theta);
	rotateMat[0][2] = vector[2]*vector[0]*(1 - cos(theta)) - vector[1]*sin(theta);	rotateMat[1][2] = vector[2]*vector[1]*(1 - cos(theta)) + vector[0]*sin(theta);	rotateMat[2][2] = pow(vector[2], 2)*(1 - cos(theta)) + cos(theta);

	multMatrix4x4(curMatrix, rotateMat);
}

/********************************************************************************
 * 									drawLine()									*
 *******************************************************************************/
//This function takes in two points, transforms them using the current transform matrix,
//	and then draws the line.
void drawLine(int x1, int y1, int z1, int x2, int y2, int z2)
{
	//printf("%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n\n", curMatrix[0][0], curMatrix[1][0], curMatrix[2][0], curMatrix[3][0], curMatrix[0][1], curMatrix[1][1], curMatrix[2][1], curMatrix[3][1], curMatrix[0][2], curMatrix[1][2], curMatrix[2][2], curMatrix[3][2], curMatrix[0][3], curMatrix[1][3], curMatrix[2][3], curMatrix[3][3]);

	myVertex v1 = {x1, y1, z1, 1};
	myVertex v2 = {x2, y2, z2, 1};


	//TRANSFORM POINTS!!!
	double* tv1 = transformVertex(curMatrix, v1);
	double* tv2 = transformVertex(curMatrix, v2);


	tv1 = transformVertex(vrcMatrix, tv1);
	tv2 = transformVertex(vrcMatrix, tv2);

	tv1 = transformVertex(projMatrix, tv1);
	tv2 = transformVertex(projMatrix, tv2);

	double* cv1 = tv1;
	double* cv2 = tv2;
	clipLine(tv1, tv2, -1, 1, -1, 1, cv1, cv2);

	tv1 = transformVertex(viewMatrix, tv1);
	tv2 = transformVertex(viewMatrix, tv2);

	/*
	int i;
	for(i = 0; i < 4; i++)
	{
		printf("%f\t%f\n", tv1[i], tv2[i]);
	}
	printf("\n");
	*/

	glBegin(GL_LINES);
		glVertex3d(tv1[0], tv1[1], tv1[2]);
		glVertex3d(tv2[0], tv2[1], tv2[2]);
	glEnd();
}

/********************************************************************************
 * 									drawLineStripv()							*
 *******************************************************************************/
//This function draws a line strip from a set of points.
void drawLineStripv(double *vector, int numVertices)
{
	int i;
	for(i = 0; i < numVertices-1; i+=3)
	{
		drawLine(vector[i], vector[i+1], vector[i+2], vector[i+3], vector[i+4], vector[i+5]);
	}
}

/********************************************************************************
 * 									drawPolyv()									*
 *******************************************************************************/
//This function draws a polygon from a set of points.
void drawPolyv(double *vector, int numVertices)
{
	int i, j;
	for(i = 0; i < numVertices-1; i++)
	{
		j = i*3;
		drawLine(vector[j], vector[j+1], vector[j+2], vector[j+3], vector[j+4], vector[j+5]);
	}

	j = i*3;
	drawLine(vector[j], vector[j+1], vector[j+2], vector[0], vector[1], vector[2]);
}


//PROJCTION
//
//
//STARTS
//
//
//HERE


/********************************************************************************
 * 									setVRC()									*
 *******************************************************************************/
//This function defines the VRC, or Viewing Reference Coordinate system. It is defined by theee
// points, VRP, At, and UP.
void setVRC(double vrpx, double vrpy, double vrpz, double atx, double aty, double atz, double upx, double upy, double upz)
{
	vrp[0] = vrpx;
	vrp[1] = vrpy;
	vrp[2] = vrpz;

	vup[0] = upx;
	vup[1] = upy;
	vup[2] = upz;

	//Calculate n vector: VRP - AT
	double nAxis[3];
	nAxis[0] = vrpx - atx;
	nAxis[1] = vrpy - aty;
	nAxis[2] = vrpz - atz;

	//Normalize
	double length = sqrt(pow(nAxis[0], 2) + pow(nAxis[1], 2) + pow(nAxis[2], 2));
	nAxis[0] /= length;
	nAxis[1] /= length;
	nAxis[2] /= length;

	//Calculate u vector: VPN X VUP
	double uAxis[3];
	uAxis[0] = nAxis[1]*upz - nAxis[2]*upy;
	uAxis[1] = nAxis[0]*upz - nAxis[2]*upx;
	uAxis[2] = nAxis[0]*upy - nAxis[1]*upx;

	//Calculate v vector: n X u
	double vAxis[3];
	vAxis[0] = nAxis[1]*uAxis[2] - nAxis[2]*uAxis[1];
	vAxis[1] = nAxis[0]*uAxis[2] - nAxis[2]*uAxis[0];
	vAxis[2] = nAxis[0]*uAxis[1] - nAxis[1]*uAxis[0];


	//Set to projection matrix
	vrcMatrix[0][0] = uAxis[0];
	vrcMatrix[0][1] = uAxis[1];
	vrcMatrix[0][2] = uAxis[2];

	vrcMatrix[1][0] = vAxis[0];
	vrcMatrix[1][1] = vAxis[1];
	vrcMatrix[1][2] = vAxis[2];

	vrcMatrix[2][0] = nAxis[0];
	vrcMatrix[2][1] = nAxis[1];
	vrcMatrix[2][2] = nAxis[2];

	//Translate VRP to origin
	myMat translate = {{1, 0, 0, -vrpx},
					   {0, 1, 0, -vrpy},
					   {0, 0, 1, -vrpz},
					   {0, 0, 0, 1}};

	multMatrix4x4(vrcMatrix, translate);
}

/********************************************************************************
 * 									setProjection()								*
 *******************************************************************************/
//This function creates the matrix that all points get multiplied by during projection.
void setProjection(int type, double prpx, double prpy, double prpz, double width, double height, double front, double back)
{
	TYPE = type;

	//PRP = COP
	//Create projection matrix using general case.
	//dx, dy, dz represent unit vector from zp (0,0, zp) to COP.
	//Q defined as distance from zp to COP.
	double dx, dy, dz, zp, Q;
	double dop[3];

	if(TYPE==PARALLEL)
	{
		zp = 0;
	}
	else
	{
		zp = (back-front/2)+front;
	}

	dx = -prpx;
	dy = -prpy;
	dz = zp-prpz;

	Q = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
	dx /= Q;
	dy /= Q;
	dz /= Q;

	dop[0] = (viewMatrix[0][0]/2) - prpx;
	dop[1] = (viewMatrix[1][1]/2) - prpy;
	dop[2] = -prpz;

	myMat proj = {{1, 0, -dx/dz, zp*(dx/dz)},
				  {0, 1, -dy/dz, zp*(dy/dz)},
				  {0, 0, -zp/(Q*dz), pow(zp, 2)/(Q*dz)+zp},
				  {0, 0, -1/(Q*dz), zp/(Q*dz)+1}};

	//Transform into Canonical View Volume
	//Rotate so n -> z, v -> y, and u -> x
	double rz[3] = {vrcMatrix[0][0], vrcMatrix[0][1], vrcMatrix[0][2]}, rx[3], ry[3];

	//rx = (vup x rz)/abs(vup x rz)
	rx[0] = vup[1]*rz[2] - vup[2]*rz[1];	//yz - zy
	rx[1] = vup[0]*rz[2] - vup[2]*rz[0];	//xz - zx
	rx[2] = vup[0]*rz[1] - vup[1]*rz[0];	//xy - yx

	double rxMag = sqrt(pow(rx[0], 2) + pow(rx[1], 2) + pow(rx[2], 2));
	rx[0] /= rxMag;
	rx[1] /= rxMag;
	rx[2] /= rxMag;

	//ry = rz x rx
	ry[0] = rz[1]*rx[2] - rz[2]*rx[1];	//yz - zy
	ry[1] = rz[0]*rx[2] - rz[2]*rx[0];	//xy - zx
	ry[2] = rz[0]*rx[1] - rz[1]*rx[0];	//xy - yx

	//Rotate so VRC aligns with global axis
	myMat rotate = {{rx[0], rx[1], rx[2], 0},
					{ry[0], ry[1], ry[2], 0},
					{rz[0], rz[1], rz[2], 0},
					{0, 0, 0, 1}};

	//Shear to bring DOP into z axis
	double hx, hy;
	hx = -dop[0]/dop[2];
	hy = -dop[1]/dop[2];

	myMat shear = {{1, 0, hx, 0},
				   {0, 1, hy, 0},
				   {0, 0, 1, 0},
				   {0, 0, 0, 1}};

	double uMax = width/2, uMin = -width/2, vMax = height/2, vMin = -height/2;
	if(TYPE == PARALLEL)
	{
		//Translate front center of View Volume to origin
		myMat translatePar = {{1, 0, 0, (-uMax+uMin)/2},
							  {0, 1, 0, (-vMax+vMin)/2},
							  {0, 0, 1, -front},
							  {0, 0, 0, 1}};

		//Scale to 2x2x1
		myMat scalePar = {{2/(uMax-uMin), 0, 0, 0},
						  {0, 2/(vMax-vMin), 0, 0},
						  {0, 0, 1/(front-back), 0},
						  {0, 0, 0, 1}};

		multMatrix4x4(projMatrix, proj);
		multMatrix4x4(projMatrix, rotate);
		multMatrix4x4(projMatrix, shear);
		multMatrix4x4(projMatrix, translatePar);
		multMatrix4x4(projMatrix, scalePar);

	}
	else
	{
		//Translate COP -> origin
		//Technically, this step goes before Shear to bring DOP into z axis
		myMat translatePer = {{1, 0, 0, -prpx},
							  {0, 1, 0, -prpy},
							  {0, 0, 1, -prpz},
							  {0, 0, 0, 1}};

		//Scale into Canonical view volume
		double vrpPrimeZ = vrcMatrix[2][2] - prpz;
		myMat scalePer = {{(2*vrpPrimeZ)/((uMax-uMin)*(vrpPrimeZ+back)), 0, 0, 0},
						  {0, (2*vrpPrimeZ)/((vMax-vMin)*(vrpPrimeZ+back)), 0, 0},
						  {0, 0, -1/(vrpPrimeZ+back), 0},
						  {0, 0, 0, 1}};

		//Change into Parallel Canonical view volume
		myMat M = {{1, 0, 0, 0},
				   {0, 1, 0, 0},
				   {0, 0, 1/(1+front), -front/(1+front)},
				   {0, 0, -1, 0}};

		multMatrix4x4(projMatrix, proj);
		multMatrix4x4(projMatrix, rotate);
		multMatrix4x4(projMatrix, translatePer);
		multMatrix4x4(projMatrix, shear);
		multMatrix4x4(projMatrix, scalePer);
		multMatrix4x4(projMatrix, M);
	}
}

/********************************************************************************
 * 									setView()									*
 *******************************************************************************/
//This function sets up the viewport that all objects will be drawn into.
void setView(int xMin, int yMin, int xMax, int yMax)
{
	//Similar to glViewport
	//Translate to xMin, yMin
	//Scale to (xMax - xMin) by (yMax - yMin)
	//x = xMin, y = yMin; width and height derived from arguments.
	double windowWidth = glutGet(GLUT_WINDOW_WIDTH), windowHeight = glutGet(GLUT_WINDOW_HEIGHT);
	double width = (xMax - xMin)/windowWidth;
	double height = (yMax - yMin)/windowHeight;

	viewMatrix[0][0] = width;
	viewMatrix[1][1] = height;
	viewMatrix[3][0] = xMin;
	viewMatrix[3][1] = yMin;
}

/********************************************************************************
 * 									clipLine()									*
 *******************************************************************************/
//This function clips a line defined by points p1 and p2 against the canonical
// view volume using the Cyrus-Beck 2D clipping algorithm.
int clipLine(myVertex p1, myVertex p2, double xMin, double xMax, double yMin, double yMax, myVertex p3, myVertex p4)
{
	//0 = Top; 1 = Left; 2 = Right; 3 = Bottom
	myVertex normals[4] = {{0, 1, 0, 0}, {-1, 0, 0, 0}, {1, 0, 0, 0}, {0, -1, 0, 0}};
	myVertex pE[2] = {{xMin, yMin, 0, 1}, {xMax, yMax, 0, 1}};

	double t[4];
	double sign[4];

	//Calculate t values of intersection of each edge
	int i;
	for(i = 0; i < 4; i++)
	{
		myVertex p0pe;
		p0pe[0] = p1[0]-pE[i%2][0];
		p0pe[1] = p1[1]-pE[i%2][1];
		p0pe[2] = p1[2]-pE[i%2][2];
		p0pe[3] = p1[3]-pE[i%2][3];

		myVertex p1p0;
		p1p0[0] = p2[0]-p1[0];
		p1p0[1] = p2[1]-p1[1];
		p1p0[2] = p2[2]-p1[2];
		p1p0[3] = p2[3]-p1[3];

		t[i] = (normals[i][0]*p0pe[0] + normals[i][1]*p0pe[1] + normals[i][2]*p0pe[2])/-(normals[i][0]*p1p0[0] + normals[i][1]*p1p0[1] + normals[i][2]*p1p0[2]);
		sign[i] = normals[i][0]*p1p0[0] + normals[i][1]*p1p0[1] + normals[i][2]*p1p0[2];
	}

	//Sort t values
	int j;
	double temp;
	for(i = 3; i > 0; i--)
	{
		for(j = 1; j <= i; j++)
		{
			if(t[j-1] > t[j])
			{
				temp = t[j-1];
				t[j-1] = t[j];
				t[j] = temp;

				temp = sign[j-1];
				sign[j-1] = sign[j];
				sign[j] = temp;
			}
		}
	}

	double tE, tL;
	//Classify points as pE or pL
	for(i = 0; i < 4; i++)
	{
		if(sign[i]<0)
		{
			if(t[i]>tE)
			{
				tE = t[i];
			}
		}
		else
		{
			if(t[i]<tL)
			{
				tL = t[i];
			}
		}
	}

	//If tE > tL, line is completely exterior to clip rectangle; trivially reject
	if(tE > tL)
	{
		return 0;
	}

	//Return the new endpoints in p3 and p4.
	if(0<=tE && tE<=1)
	{
		p3[0] = p1[0] + tE*(p2[0] - p1[0]);
		p3[1] = p1[1] + tE*(p2[1] - p1[1]);
		p3[2] = p1[2] + tE*(p2[2] - p1[2]);
		p3[3] = p1[3] + tE*(p2[3] - p1[3]);
	}
	else
	{
		p3 = p1;
	}

	//Same for second point
	if(0<=tL && tL<=1)
	{
		p4[0] = p1[0] + tL*(p2[0] - p1[0]);
		p4[1] = p1[1] + tL*(p2[1] - p1[1]);
		p4[2] = p1[2] + tL*(p2[2] - p1[2]);
		p4[3] = p1[3] + tL*(p2[3] - p1[3]);
	}
	else
	{
		p4 = p2;
	}

	return 1;
}

#endif /* MYAPI_H_ */
