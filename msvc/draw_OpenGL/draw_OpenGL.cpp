/* draw_OpenGL: postprocessing tool */

#include <windows.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Gl\glew.h>
#include <Gl\freeglut.h>
#include "GL\glui.h"
using namespace std;

/*************** Definition of global variables ***************/
INT *SATOMS, *ATOMSQC, *BONDSQC, *TRIANGLES; // data are set as global to be accessible to all OpenGL functions
FLOAT *BONDS, *R, *Z, *K, *TIME, *SIZED;
int TINDX, VARIABLESSW, RAMPID, ATOMSID, BONDSID, TRIANGLESSW, ATOMSSW, BONDSSW, MAIN_WINDOW, SS;
float MINZ, MAXZ, MINK, MAXK, DEFSCALE, ISW, SW;
float ZOOM, POSITION_X, POSITION_Y, LAST_X, LAST_Y;
GLUI *glui;

/******************** Function declarations *******************/
void colormap(float s, GLfloat* colours, int sw); // Build colormap
void setViewport(GLint left, GLint right, GLint bottom, GLint top); // Set Viewport
void myInit(void); // My Init
void myKeyboard(unsigned char theKey, int mouseX, int mouseY); // My Keyboard
void myMouse(int mouseButton, int mouseButtonState, int mouseX, int mouseY); // My Mouse
void myMotion(int mouseX, int mouseY); // My Motion
void myReshape(int w, int h); // My Reshape
void processMenuEvents(int option); // My Menu
void myDisplay(); // My Display
void setWindow(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top); // Set Window
void findMaxMin(); // Find MINZ, MAXZ, MINK, MAXK
void plotAtoms(); // Plot atoms
void plotBonds(); // Plot bonds
void plotTriangles(); // Plot triangles

/************************ Main program ************************/
void main(){

	// Create named shared memory and assign buffers
	HANDLE hsatoms, hatomsQC, hbonds, hbondsQC, htriangles, hR, hZ, hK, hTime, hSize;
	hsatoms = OpenFileMapping(FILE_MAP_READ, FALSE, "satoms");
	hatomsQC = OpenFileMapping(FILE_MAP_READ, FALSE, "atomsQC");
	hbonds = OpenFileMapping(FILE_MAP_READ, FALSE, "bonds");
	hbondsQC = OpenFileMapping(FILE_MAP_READ, FALSE, "bondsQC");
	htriangles = OpenFileMapping(FILE_MAP_READ, FALSE, "triangles");
	hR = OpenFileMapping(FILE_MAP_READ, FALSE, "R");
	hZ = OpenFileMapping(FILE_MAP_READ, FALSE, "Z");
	hK = OpenFileMapping(FILE_MAP_READ, FALSE, "K");
	hTime = OpenFileMapping(FILE_MAP_READ, FALSE, "Time");
	hSize = OpenFileMapping(FILE_MAP_READ, FALSE, "size");
	if ((hsatoms == NULL) || (hatomsQC == NULL) || (hbonds == NULL) || (hbondsQC == NULL) || (htriangles == NULL) || (hR == NULL) || (hZ == NULL) || (hK == NULL) || (hTime == NULL) || (hSize == NULL)){
		cout << "Could not create file mapping objects." << endl;
		system("pause");
		return;
	}
	SATOMS = (INT*)MapViewOfFile(hsatoms, FILE_MAP_READ, 0, 0, FALSE);
	ATOMSQC = (INT*)MapViewOfFile(hatomsQC, FILE_MAP_READ, 0, 0, FALSE);
	BONDS = (FLOAT*)MapViewOfFile(hbonds, FILE_MAP_READ, 0, 0, FALSE);
	BONDSQC = (INT*)MapViewOfFile(hbondsQC, FILE_MAP_READ, 0, 0, FALSE);
	TRIANGLES = (INT*)MapViewOfFile(htriangles, FILE_MAP_READ, 0, 0, FALSE);
	R = (FLOAT*)MapViewOfFile(hR, FILE_MAP_READ, 0, 0, FALSE);
	Z = (FLOAT*)MapViewOfFile(hZ, FILE_MAP_READ, 0, 0, FALSE);
	K = (FLOAT*)MapViewOfFile(hK, FILE_MAP_READ, 0, 0, FALSE);
	TIME = (FLOAT*)MapViewOfFile(hTime, FILE_MAP_READ, 0, 0, FALSE);
	SIZED = (FLOAT*)MapViewOfFile(hSize, FILE_MAP_READ, 0, 0, FALSE);
	if ((SATOMS == NULL) || (ATOMSQC == NULL) || (BONDS == NULL) || (BONDSQC == NULL) || (TRIANGLES == NULL) || (R == NULL) || (Z == NULL) || (K == NULL) || (TIME == NULL) || (SIZED == NULL)){
		cout << "Could not map view of files." << endl;
		CloseHandle(hsatoms);
		CloseHandle(hatomsQC);
		CloseHandle(hbonds);
		CloseHandle(hbondsQC);
		CloseHandle(htriangles);
		CloseHandle(hR);
		CloseHandle(hZ);
		CloseHandle(hK);
		CloseHandle(hTime);
		CloseHandle(hSize);
		system("pause");
		return;
	}

	// Initialize MINZ, MAXZ, MINK, MAXK
	findMaxMin();

	// Set up size
	ISW = (float) 1.3*max(SIZED[2], SIZED[3]);

	// Initialize GLUT
	int ac = 1;
	char* av = "MATLAB test";
	glutInit(&ac, &av);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

	// Get screen size
	SS = min((int) GetSystemMetrics(SM_CXSCREEN), GetSystemMetrics(SM_CYSCREEN));
	glutInitWindowSize(SS, SS);
	glutInitWindowPosition(0, 0);
	MAIN_WINDOW = glutCreateWindow("Draw simulation results");
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);
	glutMotionFunc(myMotion);
	glutMouseFunc(myMouse);
	glutKeyboardFunc(myKeyboard);
	myInit();
	SW = ISW;
	setWindow(-ISW, ISW, -ISW, ISW);

	// Set GLUI window
	glui = GLUI_Master.create_glui("GLUI", 0, SS, 0);
	glui->add_statictext("Plot options");
	glui->add_separator();

	// Draw atoms
	GLUI_Checkbox* checkbox_atoms = glui->add_checkbox("Draw atoms", &ATOMSSW);

	// Draw bonds
	GLUI_Checkbox* checkbox_bonds = glui->add_checkbox("Draw bonds", &BONDSSW);
	checkbox_bonds->set_int_val(1);

	// Draw triangles
	GLUI_Checkbox* checkbox_triangles = glui->add_checkbox("Draw triangles", &TRIANGLESSW);
	glui->add_separator();

	// Choose atoms to draw
	GLUI_Listbox* listbox_atoms = glui->add_listbox("Atoms     ", &ATOMSID, NULL);
	listbox_atoms->add_item(1, "All atoms");
	listbox_atoms->add_item(2, "Satoms");
	listbox_atoms->add_item(3, "Repatoms");
	listbox_atoms->set_int_val(1);

	// Choose variable
	GLUI_Listbox* listbox_variable = glui->add_listbox("Variable  ", &VARIABLESSW, NULL);
	listbox_variable->add_item(1, "Plastic deformation");
	listbox_variable->add_item(2, "Cum. plast. def.");
	listbox_variable->set_int_val(1);

	// Choose color ramp
	GLUI_Listbox* listbox_color = glui->add_listbox("Clr ramp  ", &RAMPID, NULL);
	listbox_color->add_item(1, "RGB bipolar");
	listbox_color->add_item(2, "Thermal");
	listbox_color->add_item(3, "JET");
	listbox_color->add_item(4, "RWB");
	listbox_color->set_int_val(3);

	// Choose bonds to draw
	GLUI_Listbox* listbox_bonds = glui->add_listbox("Bonds     ", &BONDSID, NULL);
	listbox_bonds->add_item(1, "All bonds");
	listbox_bonds->add_item(2, "Sbonds");
	listbox_bonds->set_int_val(1);
	glui->add_separator();

	// Time instant
	GLUI_Spinner* spinner_time = glui->add_spinner("Time step", GLUI_SPINNER_INT, &TINDX);
	spinner_time->set_int_limits(0, (int)(max(TIME[0], TIME[1]) - 1));

	// Scale
	GLUI_Spinner* spinner_scale = glui->add_spinner("Scale     ", GLUI_SPINNER_FLOAT, &DEFSCALE);
	spinner_scale->set_int_limits(0, 100);
	glui->add_separator();

	// Translation Z
	GLUI_Translation* Translation_z = new GLUI_Translation(glui, "Zoom", GLUI_TRANSLATION_Z, &ZOOM);
	Translation_z->set_speed((GLfloat)(5*ISW / SS)); // invariant z-scaling speed

	// Link glut to glui
	glui->set_main_gfx_window(MAIN_WINDOW);

	// Launch the main loop
	glutMainLoop();

	// Unmap buffers
	UnmapViewOfFile(SATOMS);
	UnmapViewOfFile(ATOMSQC);
	UnmapViewOfFile(BONDS);
	UnmapViewOfFile(BONDSQC);
	UnmapViewOfFile(TRIANGLES);
	UnmapViewOfFile(R);
	UnmapViewOfFile(Z);
	UnmapViewOfFile(K);
	UnmapViewOfFile(TIME);
	UnmapViewOfFile(SIZED);
	CloseHandle(hsatoms);
	CloseHandle(hatomsQC);
	CloseHandle(hbonds);
	CloseHandle(hbondsQC);
	CloseHandle(htriangles);
	CloseHandle(hR);
	CloseHandle(hZ);
	CloseHandle(hK);
	CloseHandle(hTime);
	CloseHandle(hSize);

	return;
}

/******************** Function definitions ********************/

/*<<<<<<<<<<<<<<<<<<< Set Window >>>>>>>>>>>>>>>>>>>*/
void setWindow(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top){
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(left, right, bottom, top);
}

/*<<<<<<<<<<<<<<<<<<< Set Viewport >>>>>>>>>>>>>>>>>>>*/
void setViewport(GLint left, GLint right, GLint bottom, GLint top){
	glViewport(left, bottom, right - left, top - bottom);
}

/*<<<<<<<<<<<<<<<<<<< My Init >>>>>>>>>>>>>>>>>>>*/
void myInit(void){
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_NOTEQUAL, 0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(6.0);
	glLineWidth(1.0);
}

/*<<<<<<<<<<<<<<<<<<< My Keyboard >>>>>>>>>>>>>>>>>>>*/
void myKeyboard(unsigned char theKey, int mouseX, int mouseY){
	int NT = (int)max(TIME[0], TIME[1]);
	switch (theKey){
	case '+':
		if (TINDX < (NT - 1)){
			TINDX++;
		}
		break;
	case '-':
		if (TINDX > 0){
			TINDX--;
		}
		break;
	default:
		break;
	}
	glui->sync_live();
	glutPostRedisplay();
}

/*<<<<<<<<<<<<<<<<<<< My Mouse >>>>>>>>>>>>>>>>>>>*/
void myMouse(int mouseButton, int mouseButtonState, int mouseX, int mouseY){
	if ((mouseButton == GLUT_LEFT_BUTTON) && (mouseButtonState == GLUT_DOWN)){
		LAST_X = (float)mouseX;
		LAST_Y = (float)mouseY;
	}
}

/*<<<<<<<<<<<<<<<<<<< My Reshape >>>>>>>>>>>>>>>>>>>*/
void myReshape(int w, int h){
	SS = min(w,h);
	glutReshapeWindow(SS, SS);
	glViewport(0, 0, SS, SS);	
}

/*<<<<<<<<<<<<<<<<<<< My Motion >>>>>>>>>>>>>>>>>>>*/
void myMotion(int mouseX, int mouseY){
	POSITION_X += (float)((mouseX - LAST_X) * 2 * SW / SS);
	POSITION_Y += (float)((-mouseY + LAST_Y) * 2 * SW / SS);
	LAST_X = (float)mouseX;
	LAST_Y = (float)mouseY;
	glutPostRedisplay();
}

/*<<<<<<<<<<<<<<<<<<< My Menu >>>>>>>>>>>>>>>>>>>*/
void processMenuEvents(int option){
	switch (option){
	case 1:
		VARIABLESSW = 1;
		break;
	case 2:
		VARIABLESSW = 2;
		break;
	}
	glutPostRedisplay();
}

/*<<<<<<<<<<<<<<<<<<< My Display >>>>>>>>>>>>>>>>>>>*/
void myDisplay(){ // plot the structure
	SW = (ISW - ZOOM);
	SW = (ISW - ZOOM);
	setWindow(-SW, SW, -SW, SW);
	glTranslatef((GLfloat)POSITION_X, (GLfloat)POSITION_Y, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// Plot triangles
	if (TRIANGLESSW == 1){
		plotTriangles();
	}

	// Plot atoms
	if (ATOMSSW == 1){
		plotAtoms();
	}
	
	// Plot bonds
	if (BONDSSW == 1){
		plotBonds();
	}
	
	// Draw
	glutSwapBuffers();
}

/*<<<<<<<<<<<<<<<<<<< Plot atoms >>>>>>>>>>>>>>>>>>>*/
void plotAtoms(){
	int i, alpha;
	int MR = (int)R[0];
	int Na = (int)R[0] / 2;
	int NQC = (int)ATOMSQC[0];
	int Ns = (int)SATOMS[0];
	switch (ATOMSID){
	case 1:
		glBegin(GL_POINTS);
		for (i = 0; i < Na; i++){
			glColor4f(0.6f, 0.6f, 0.6f, 0.75f);
			glVertex2f((GLfloat)(R[2 * i + 2] + DEFSCALE*(R[2 * i + TINDX*MR + 2] - R[2 * i + 2])), (GLfloat)(R[2 * i + 1 + 2] + DEFSCALE*(R[2 * i + 1 + TINDX*MR + 2] - R[2 * i + 1 + 2])));
		}
		glEnd();
		return;
	case 2:
		glBegin(GL_POINTS);
		for (i = 0; i < Ns; i++){
			glColor4f(0.6f, 0.6f, 0.6f, 0.75f);
			alpha = (int)SATOMS[i + 2] - 1; // MATLAB indexing
			glVertex2f((GLfloat)(R[2 * alpha + 2] + DEFSCALE*(R[2 * alpha + TINDX*MR + 2] - R[2 * alpha + 2])), (GLfloat)(R[2 * alpha + 1 + 2] + DEFSCALE*(R[2 * alpha + 1 + TINDX*MR + 2] - R[2 * alpha + 1 + 2])));
		}
		glEnd();
		return;
	case 3:
		glBegin(GL_POINTS);
		for (i = 0; i < NQC; i++){
			glColor4f(0.6f, 0.6f, 0.6f, 0.75f);
			alpha = (int)ATOMSQC[i + 2] - 1; // MATLAB indexing
			glVertex2f((GLfloat)(R[2 * alpha + 2] + DEFSCALE*(R[2 * alpha + TINDX*MR + 2] - R[2 * alpha + 2])), (GLfloat)(R[2 * alpha + 1 + 2] + DEFSCALE*(R[2 * alpha + 1 + TINDX*MR + 2] - R[2 * alpha + 1 + 2])));
		}
		glEnd();
		return;
	}
}

/*<<<<<<<<<<<<<<<<<<< Plot bonds >>>>>>>>>>>>>>>>>>>*/
void plotBonds(){
	int MR = (int)R[0];

	// Plot colormap
	glBegin(GL_POINTS);
	GLfloat colours[3];
	float param;
	for (param = 0; param <= 1; param += 0.002f){
		colormap(param, colours, RAMPID);
		glColor4f(colours[0], colours[1], colours[2], 0.75f);
		glVertex2f((GLfloat)(-ISW*0.9), (GLfloat)(0.9*ISW*(2 * param - 1)));
	}
	glEnd();

	// Plot bonds in deformed configuration
	float min, max;
	FLOAT *prD;
	int MD;
	switch (VARIABLESSW){ // find MIN, MAX
	case 1:
		prD = Z;
		MD = (int)Z[0];
		min = MINZ;
		max = MAXZ;
		break;
	case 2:
		prD = K;
		MD = (int)K[0];
		min = MINK;
		max = MAXK;
		break;
	}
	int alpha, beta, bond, i;
	int MB = (int)BONDS[0];
	int MBQC = (int)BONDSQC[0];
	switch (BONDSID){
	case 1:
		glBegin(GL_LINES);
		for (i = 0; i < MB; i++){
			alpha = (int)BONDS[i + MB * 0 + 2] - 1; // atom alpha, MATLAB indexing
			beta = (int)BONDS[i + MB * 1 + 2] - 1; // atom beta, MATLAB indexing

			// Set color according to value param
			param = (prD[i + TINDX*MD + 2] - min) / (max - min);
			colormap(param, colours, RAMPID);
			glColor4f(colours[0], colours[1], colours[2], 0.75f);
			glVertex2f((GLfloat)(R[2 * alpha + 2] + DEFSCALE*(R[2 * alpha + TINDX*MR + 2] - R[2 * alpha + 2])), (GLfloat)(R[2 * alpha + 1 + 2] + DEFSCALE*(R[2 * alpha + 1 + TINDX*MR + 2] - R[2 * alpha + 1 + 2])));
			glVertex2f((GLfloat)(R[2 * beta + 2] + DEFSCALE*(R[2 * beta + TINDX*MR + 2] - R[2 * beta + 2])), (GLfloat)(R[2 * beta + 1 + 2] + DEFSCALE*(R[2 * beta + 1 + TINDX*MR + 2] - R[2 * beta + 1 + 2])));
		}
		glEnd();
		break;
	case 2:
		glBegin(GL_LINES);
		for (i = 0; i < MBQC; i++){
			bond = (int)BONDSQC[i + 2] - 1; // repbond, MATLAB indexing
			alpha = (int)BONDS[bond + MB * 0 + 2] - 1; // atom alpha, MATLAB indexing
			beta = (int)BONDS[bond + MB * 1 + 2] - 1; // atom beta, MATLAB indexing

			// Set color according to value param
			param = (prD[bond + TINDX*MD + 2] - min) / (max - min);
			colormap(param, colours, RAMPID);
			glColor4f(colours[0], colours[1], colours[2], 0.75f);
			glVertex2f((GLfloat)(R[2 * alpha + 2] + DEFSCALE*(R[2 * alpha + TINDX*MR + 2] - R[2 * alpha + 2])), (GLfloat)(R[2 * alpha + 1 + 2] + DEFSCALE*(R[2 * alpha + 1 + TINDX*MR + 2] - R[2 * alpha + 1 + 2])));
			glVertex2f((GLfloat)(R[2 * beta + 2] + DEFSCALE*(R[2 * beta + TINDX*MR + 2] - R[2 * beta + 2])), (GLfloat)(R[2 * beta + 1 + 2] + DEFSCALE*(R[2 * beta + 1 + TINDX*MR + 2] - R[2 * beta + 1 + 2])));
		}
		glEnd();
		break;
	}

	// Render text
	glColor4f(0.6f, 0.6f, 0.6f, 0.75f);

	// Min
	stringstream temptextmin;
	string textmin;
	temptextmin << min;
	textmin = temptextmin.str();
	glRasterPos2f((GLfloat)(-ISW*0.9), (GLfloat)(-0.9*ISW));
	for (i = 0; i < (int)textmin.length(); i++)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (const char)textmin[i]);
	
	// Max
	stringstream temptextmax;
	string textmax;
	temptextmax << max;
	textmax = temptextmax.str();
	glRasterPos2f((GLfloat)(-ISW*0.9), (GLfloat)(0.9*ISW));
	for (i = 0; i < (int)textmax.length(); i++)
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (const char)textmax[i]);
}

/*<<<<<<<<<<<<<<<<<<< Plot triangles >>>>>>>>>>>>>>>>>>>*/
void plotTriangles(){
	int i, alpha, beta, gamma;
	int MR = (int)R[0];
	int Nt = (int)TRIANGLES[0];

	// Plot triangles in deformed configuration
	for (i = 0; i < Nt; i++){
		glBegin(GL_LINE_LOOP);
		alpha = (int)TRIANGLES[i + Nt * 0 + 2] - 1; // atom alpha, MATLAB indexing
		beta = (int)TRIANGLES[i + Nt * 1 + 2] - 1; // atom beta, MATLAB indexing
		gamma = (int)TRIANGLES[i + Nt * 2 + 2] - 1; // atom beta, MATLAB indexing
		glColor4f(0.6f, 0.6f, 0.6f, 0.75f);
		glVertex2f((GLfloat)(R[2 * alpha + 2] + DEFSCALE*(R[2 * alpha + TINDX*MR + 2] - R[2 * alpha + 2])), (GLfloat)(R[2 * alpha + 1 + 2] + DEFSCALE*(R[2 * alpha + 1 + TINDX*MR + 2] - R[2 * alpha + 1 + 2])));
		glVertex2f((GLfloat)(R[2 * beta + 2] + DEFSCALE*(R[2 * beta + TINDX*MR + 2] - R[2 * beta + 2])), (GLfloat)(R[2 * beta + 1 + 2] + DEFSCALE*(R[2 * beta + 1 + TINDX*MR + 2] - R[2 * beta + 1 + 2])));
		glVertex2f((GLfloat)(R[2 * gamma + 2] + DEFSCALE*(R[2 * gamma + TINDX*MR + 2] - R[2 * gamma + 2])), (GLfloat)(R[2 * gamma + 1 + 2] + DEFSCALE*(R[2 * gamma + 1 + TINDX*MR + 2] - R[2 * gamma + 1 + 2])));
		glEnd();
	}
}

/*<<<<<<<<<<<<<<<<<<< Build colormap >>>>>>>>>>>>>>>>>>>*/
void colormap(float s, GLfloat* colours, int sw){
	switch (sw){
	case 1:
		if (s <= 0.25){
			colours[0] = 0;
			colours[1] = (GLfloat)(-s / 0.25 + 1);
			colours[2] = 1;
		}
		else if ((s > 0.25) && (s <= 0.5)){
			colours[0] = 0;
			colours[1] = 0;
			colours[2] = (GLfloat)(-s / 0.25 + 2);
		}
		else if ((s > 0.5) && (s <= 0.75)){
			colours[0] = (GLfloat)(s / 0.25 - 2);
			colours[1] = 0;
			colours[2] = 0;
		}
		else if ((s > 0.75) && (s <= 1)){
			colours[0] = 1;
			colours[1] = (GLfloat)(s / 0.25 - 3);
			colours[2] = 0;
		}
		break;
	case 2:
		if (s <= 0.333333333333333){
			colours[0] = (GLfloat)(3 * s);
			colours[1] = 0;
			colours[2] = 0;
		}
		else if ((s > 0.333333333333333) && (s <= 0.666666666666667)){
			colours[0] = 1;
			colours[1] = (GLfloat)(3 * s - 1);
			colours[2] = 0;
		}
		else if ((s > 0.666666666666667) && (s <= 1)){
			colours[0] = 1;
			colours[1] = 1;
			colours[2] = (GLfloat)(3 * s - 2);
		}
		break;
	case 3:
		if (s <= 0.25){
			colours[0] = 0;
			colours[1] = (GLfloat)(s / 0.25);
			colours[2] = 1;
		}
		else if ((s > 0.25) && (s <= 0.5)){
			colours[0] = 0;
			colours[1] = 1;
			colours[2] = (GLfloat)(-s / 0.25 + 2);
		}
		else if ((s > 0.5) && (s <= 0.75)){
			colours[0] = (GLfloat)(s / 0.25 - 2);
			colours[1] = 1;
			colours[2] = 0;
		}
		else if ((s > 0.75) && (s <= 1)){
			colours[0] = 1;
			colours[1] = (GLfloat)(-s / 0.25 + 4);
			colours[2] = 0;
		}
		break;
	case 4:
		if (s <= 0.5){
			colours[0] = (GLfloat)(2 * s);
			colours[1] = (GLfloat)(2 * s);
			colours[2] = 1;
		}
		else if ((s > 0.5) && (s <= 1)){
			colours[0] = 1;
			colours[1] = (GLfloat)(-2 * s + 2);
			colours[2] = (GLfloat)(-2 * s + 2);
		}
		break;
	}
}

/*<<<<<<<<<<<<<<<<<<< Find MINZ, MAXZ, MINK, MAXK >>>>>>>>>>>>>>>>>>>*/
void findMaxMin(){
	int MZ = (int)Z[0];
	int NZ = (int)Z[1];
	MINZ = 0;
	MAXZ = 0;
	for (int j = 0; j < NZ; j++){
		for (int i = 0; i < MZ; i++){
			if (Z[(i + j*MZ) + 2] < MINZ){
				MINZ = Z[(i + j*MZ) + 2];
			}
			if (Z[(i + j*MZ) + 2] > MAXZ){
				MAXZ = Z[(i + j*MZ) + 2];
			}
		}
	}
	int MK = (int)K[0];
	int NK = (int)K[1];
	MINK = 0;
	MAXK = 0;
	for (int j = 0; j < NK; j++){
		for (int i = 0; i < MK; i++){
			if (K[(i + j*MK) + 2] < MINK){
				MINK = K[(i + j*MK) + 2];
			}
			if (K[(i + j*MK) + 2] > MAXK){
				MAXK = K[(i + j*MK) + 2];
			}
		}
	}
}
