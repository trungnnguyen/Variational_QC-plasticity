// refine_mesh: for given t, p, and tr_to_refine provides refined t1, p1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "mex.h"
#include "matrix.h"

/******************* Definition of constants ******************/
const int NLEPP = 100; // average length is about 7
const int MULT = 10; // memory buffer for t, p, Marker, new lengths are Nt+MULT*Nt

/******************* Definition of global variables ******************/
double TOL; // a positive number < TOL is considered as zero

/******************** Function declarations *******************/
int longest_edge(double* prp, double* prt, int id); // function returning id of the longest edge in id-th triangle
int longest_edge_neighbour(int** t, int* mark, int idLE, int Nt, int id); // function returning id of id-th triangle longest edge neighbour
void evaluate_LEPP(int**t, int Nt, int* mark, int K, int* LEPP); // function evaluating LEPP
void subdivide_triangles(int** t, int* Nt, double** p, int* Np, int* mark,
	int* LEPP); // function subdividing the last one or the two last triangles in LEPP
void correct_LEPP(int** t, int Nt, int* mark, int* LEPP); // adjust LEPP from previous step, backtrack two triangles and reevaluate LEPP tail on refined mesh

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for five input arguments
	if (nrhs != 5){
		mexErrMsgTxt("Five input arguments required.");
	}

	// Copy input variables
	double* prp = mxGetPr(prhs[0]);
	double* prt = mxGetPr(prhs[1]);
	double* prtr_to_refine = mxGetPr(prhs[2]);
	double* prSW = mxGetPr(prhs[3]);
	int Np = (int)mxGetN(prhs[0]); // number of points in p
	int Mt = (int)mxGetM(prhs[1]); // number of rows in t
	int Nt = (int)mxGetN(prhs[1]); // number of input triangles
	int Nt_refine = (int)std::max(mxGetM(prhs[2]), mxGetN(prhs[2])); // number of triangles marked for refinement
	int SW = (int)prSW[0]; // switch for reshuffling nodes
	double *prTOL = mxGetPr(prhs[4]);
	TOL = prTOL[0];





	/*<<<<<<<<<<<<<<<<<<<<<<<< Prepare working data >>>>>>>>>>>>>>>>>>>>>>>>*/

	// Allocate workt matrices out_t and out_p
	int Nout_t = Nt;
	int Nout_p = Np;
	int Npmax = Np + MULT*Np;
	int Ntmax = Nt + MULT*Nt;
	int** out_t = new int*[4];
	double** out_p = new double*[2];
	for (int i = 0; i < 4; i++){
		out_t[i] = new int[Ntmax];
	}
	out_p[0] = new double[Npmax];
	out_p[1] = new double[Npmax];

	// Allocate and assign markers from matrix prtr_to_refine
	int* Marker = new int[Ntmax];
	for (int i = 0; i < Ntmax; i++){
		Marker[i] = 0;
	}
	for (int i = 0; i < Nt_refine; i++){
		if ((int)prtr_to_refine[i] >= 1){
			Marker[(int)prtr_to_refine[i] - 1] = 1; // MATLAB indexing
		}
		else{
			mexErrMsgTxt("IDs of triangles to refine must be positive integers.");
		}
	}

	// Reshuffle indices: in t, the first two indices have to correspond to the longest edge in the triangle. Usually needs to be done only once for the initial mesh, provided arrays are consistent with this requirement
	if (SW == 1){
		for (int i = 0; i < Nt; i++){

			// Find triangle's LE
			int idLE = longest_edge(prp, prt, i);

			// Update data in out_t
			if (idLE == 2){
				out_t[0][i] = (int)prt[1 + 3 * i];
				out_t[1][i] = (int)prt[2 + 3 * i];
				out_t[2][i] = (int)prt[0 + 3 * i];
				out_t[3][i] = 0;
			}
			else if (idLE == 3){
				out_t[0][i] = (int)prt[2 + 3 * i];
				out_t[1][i] = (int)prt[0 + 3 * i];
				out_t[2][i] = (int)prt[1 + 3 * i];
				out_t[3][i] = 0;
			}
			else{
				out_t[0][i] = (int)prt[0 + 3 * i];
				out_t[1][i] = (int)prt[1 + 3 * i];
				out_t[2][i] = (int)prt[2 + 3 * i];
				out_t[3][i] = 0;
			}
		}
	}
	else{
		for (int i = 0; i < Nt; i++){

			// Update data in out_t
			out_t[0][i] = (int)prt[0 + 3 * i];
			out_t[1][i] = (int)prt[1 + 3 * i];
			out_t[2][i] = (int)prt[2 + 3 * i];
			out_t[3][i] = 0;
		}
	}
	for (int i = 0; i < Np; i++){

		// Update data in out_p
		out_p[0][i] = prp[0 + 2 * i];
		out_p[1][i] = prp[1 + 2 * i];
	}





	/*<<<<<<<<<<<<<<<<<<<<<<<< Send the data straight back to MATLAB if no refinement occurs >>>>>>>>>>>>>>>>>>>>>>>>*/
	if (Nt_refine == 0){

		// Initialize and populate
		nlhs = 2;
		plhs[0] = mxCreateDoubleMatrix(2, Np, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(4, Nt, mxREAL);
		double* prp1 = mxGetPr(plhs[0]);
		double* prt1 = mxGetPr(plhs[1]);

		// Populate p1 with data
		for (int i = 0; i < Np; i++){
			prp1[0 + 2 * i] = out_p[0][i];
			prp1[1 + 2 * i] = out_p[1][i];
		}

		// Populate t1 with data
		for (int i = 0; i < Nt; i++){
			prt1[0 + 4 * i] = out_t[0][i];
			prt1[1 + 4 * i] = out_t[1][i];
			prt1[2 + 4 * i] = out_t[2][i];
			prt1[3 + 4 * i] = -(i + 1);
		}

		// Delete temporary arrays
		for (int i = 0; i < 4; i++){
			delete[] out_t[i];
		}
		delete[] out_t;
		delete[] out_p[0];
		delete[] out_p[1];
		delete[] out_p;
		delete[] Marker;

		return;
	}





	/*<<<<<<<<<<<<<<<<<<<<<<<< Refine the mesh >>>>>>>>>>>>>>>>>>>>>>>>*/

	// Loop over all triangles marked for refinement
	int LEPP[NLEPP];
	for (int i = 0; i < Nt_refine; i++){
		int K = (int)(prtr_to_refine[i] - 1); // id of the triangle to be refined, MATLAB indexing

		// While triangle K is not refined itself and marked for deletion do
		int iter_count = 0; // LEPP iteration counter
		while (Marker[K] != -1){

			// First evaluation of LEPP
			if (iter_count == 0){
				evaluate_LEPP(out_t, Nout_t, Marker, K, LEPP); // LEPP[0] - length of the LEPP itself, LEPP[1:1+LEPP[0]] - ids of LEPP members, LEPP[LEPP[0]+2] - end switch
			}
			else{
				correct_LEPP(out_t, Nout_t, Marker, LEPP); // adjust LEPP from the previous step, backtrack one or two triangles and reevaluate LEPP tail on refined mesh
			}

			// Subdivide triangles - take the last one or two triangles in LEPP and subdivide them
			subdivide_triangles(out_t, &Nout_t, out_p, &Nout_p, Marker,
				LEPP);

			// Update counter
			iter_count++;
		}

		// Test buffer fullness
		if ((Nout_t > 0.75*Ntmax) || (Nout_p > 0.75*Npmax)){
			mexErrMsgTxt("Buffer is full, enlarge MULT and recompile.");
		}
	}





	/*<<<<<<<<<<<<<<<<<<<<<<<< Send out the data back to MATLAB >>>>>>>>>>>>>>>>>>>>>>>>*/

	// Initialize, populate, and collect the garbage
	nlhs = 2;
	plhs[0] = mxCreateDoubleMatrix(2, Nout_p, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(4, Nout_t, mxREAL);
	double* prp1 = mxGetPr(plhs[0]);
	double* prt1 = mxGetPr(plhs[1]);

	// Populate p1 with data
	for (int i = 0; i < Nout_p; i++){
		prp1[0 + 2 * i] = out_p[0][i];
		prp1[1 + 2 * i] = out_p[1][i];
	}

	// Populate t1 with data
	int count_t = 0;
	for (int i = 0; i < Nout_t; i++){
		if (Marker[i] != -1){
			prt1[0 + 4 * count_t] = out_t[0][i];
			prt1[1 + 4 * count_t] = out_t[1][i];
			prt1[2 + 4 * count_t] = out_t[2][i];
			if (out_t[3][i] == 0){
				prt1[3 + 4 * count_t] = -(i + 1);
			}
			else{
				prt1[3 + 4 * count_t] = out_t[3][i];
			}
			count_t++;
		}
	}

	// Resize allocated outputs
	mxSetN(plhs[1], count_t);

	// Delete temporary arrays
	for (int i = 0; i < 4; i++){
		delete[] out_t[i];
	}
	delete[] out_t;
	delete[] out_p[0];
	delete[] out_p[1];
	delete[] out_p;
	delete[] Marker;
}

/******************** Function definitions ********************/
// Update LEPP tail
void correct_LEPP(int** t, int Nt, int* mark, int* LEPP){
	int next;
	int end = 0;
	int count = 0;
	if (LEPP[LEPP[0] + 1] == 0){ // LEPP ends on physical boundary, backtrack one triangle
		count = LEPP[0] - 2;
	}
	else if (LEPP[LEPP[0] + 1] == 1){ // LEPP ends by two triangles, backtrack two triangles
		count = LEPP[0] - 3;
	}

	while (end == 0){
		count++;
		next = longest_edge_neighbour(t, mark, 1, Nt, LEPP[count]);
		if (next == -1){
			LEPP[count + 1] = 0; // LEPP ends on the physical boundary
			end = 1;
		}
		else if (next == LEPP[count - 1] && count > 1){
			LEPP[count + 1] = 1; // LEPP ends by two triangles
			end = 1;
		}
		else{
			LEPP[count + 1] = next;
		}
	}
	LEPP[0] = count;
}

// Subdivide the last one or two triangles in LEPP
void subdivide_triangles(int** t, int* Nt, double** p, int* Np, int* mark,
	int* LEPP){
	int idT1, idT2, id1, id2, id3, id4;

	// Decide how many triangles will be subdivided
	if (LEPP[LEPP[0] + 1] == 0){ // LEPP ends on physical boundary, subdivide one triangle

		// Bisect the longest edge of the last triangle and create new two triangles
		idT1 = LEPP[LEPP[0] + 1 - 1]; // id of the triangle to be bisected
		id1 = t[0][idT1]; // id of the longest-edge node in idT1 triangle, MATLAB indexing
		id2 = t[1][idT1]; // id of the longest-edge node in idT1 triangle, MATLAB indexing
		id3 = t[2][idT1]; // id of the remaining node in idT1 triangle, MATLAB indexing

		// Add a new node to p
		p[0][*Np] = 0.5*(p[0][id1 - 1] + p[0][id2 - 1]);
		p[1][*Np] = 0.5*(p[1][id1 - 1] + p[1][id2 - 1]);

		// Add two new triangles
		t[0][*Nt] = id3;
		t[1][*Nt] = id1;
		t[2][*Nt] = *Np + 1;

		t[0][*Nt + 1] = id2;
		t[1][*Nt + 1] = id3;
		t[2][*Nt + 1] = *Np + 1;

		// Assign parent element
		if (t[3][idT1] == 0){
			t[3][*Nt] = idT1 + 1;
			t[3][*Nt + 1] = idT1 + 1;
		}
		else{
			t[3][*Nt] = t[3][idT1];
			t[3][*Nt + 1] = t[3][idT1];
		}

		// Mark refined triangle to be deleted
		mark[*Nt] = 0; // new triangles marked as standard (no refinement, no deletion)
		mark[*Nt + 1] = 0;
		mark[idT1] = -1; // mark refined triangle idT1 for deletion
		*Nt = *Nt + 2; // two triangles added
		*Np = *Np + 1; //  one point added
	}
	else if (LEPP[LEPP[0] + 1] == 1){ // LEPP ends by two triangles, subdivide them

		// Bisect the longest edge of the last two triangles and create four new triangles
		idT1 = LEPP[LEPP[0] + 1 - 1]; // id of the first triangle to be bisected
		idT2 = LEPP[LEPP[0] + 1 - 2]; // id of the second triangle to be bisected
		id1 = t[0][idT1]; // id of the longest-edge node in idT1 triangle, MATLAB indexing
		id2 = t[1][idT1]; // id of the longest-edge node in idT1 triangle, MATLAB indexing
		id3 = t[2][idT1]; // id of the remaining node in idT1 triangle, MATLAB indexing
		id4 = t[2][idT2]; // id of the remaining node in idT2 triangle, MATLAB indexing

		// Add a new node to p
		p[0][*Np] = 0.5*(p[0][id1 - 1] + p[0][id2 - 1]);
		p[1][*Np] = 0.5*(p[1][id1 - 1] + p[1][id2 - 1]);

		// Add four new triangles
		t[0][*Nt] = id3;
		t[1][*Nt] = id1;
		t[2][*Nt] = *Np + 1;

		t[0][*Nt + 1] = id2;
		t[1][*Nt + 1] = id3;
		t[2][*Nt + 1] = *Np + 1;

		// Assign parent element
		if (t[3][idT1] == 0){
			t[3][*Nt] = idT1 + 1;
			t[3][*Nt + 1] = idT1 + 1;
		}
		else{
			t[3][*Nt] = t[3][idT1];
			t[3][*Nt + 1] = t[3][idT1];
		}

		t[0][*Nt + 2] = id1;
		t[1][*Nt + 2] = id4;
		t[2][*Nt + 2] = *Np + 1;

		t[0][*Nt + 3] = id4;
		t[1][*Nt + 3] = id2;
		t[2][*Nt + 3] = *Np + 1;

		// Assign parent element
		if (t[3][idT2] == 0){
			t[3][*Nt + 2] = idT2 + 1;
			t[3][*Nt + 3] = idT2 + 1;
		}
		else{
			t[3][*Nt + 2] = t[3][idT2];
			t[3][*Nt + 3] = t[3][idT2];
		}

		// Mark refined triangles to be deleted
		mark[*Nt] = 0; // new triangles marked as standard (no refinement, no deletion)
		mark[*Nt + 1] = 0;
		mark[*Nt + 2] = 0;
		mark[*Nt + 3] = 0;
		mark[idT1] = -1; // mark idT1 triangle for deletion
		mark[idT2] = -1; // mark idT2 triangle for deletion
		*Nt = *Nt + 4; // four triangles added
		*Np = *Np + 1; // one point added
	}
}

// Find the longest edge in an id-th triangle
int longest_edge(double* prp, double* prt, int id){
	double P1[2], P2[2], P3[2];
	double e1, e2, e3;
	int idLE = -1;

	// Triangle's vertices
	P1[0] = prp[(int)(0 + 2 * (prt[0 + 3 * id] - 1))]; // MATLAB indexing
	P1[1] = prp[(int)(1 + 2 * (prt[0 + 3 * id] - 1))];
	P2[0] = prp[(int)(0 + 2 * (prt[1 + 3 * id] - 1))];
	P2[1] = prp[(int)(1 + 2 * (prt[1 + 3 * id] - 1))];
	P3[0] = prp[(int)(0 + 2 * (prt[2 + 3 * id] - 1))];
	P3[1] = prp[(int)(1 + 2 * (prt[2 + 3 * id] - 1))];

	// Lengths of triangle's edges
	e1 = sqrt(pow(P1[0] - P2[0], 2) + pow(P1[1] - P2[1], 2));
	e2 = sqrt(pow(P2[0] - P3[0], 2) + pow(P2[1] - P3[1], 2));
	e3 = sqrt(pow(P3[0] - P1[0], 2) + pow(P3[1] - P1[1], 2));

	// Find the LE
	if ((fabs(e1 - e2) < TOL) && (fabs(e1 - e3) < TOL) &&
		(fabs(e2 - e3) < TOL)){ // all three edges are of the same length

		// Try to elect the horizontal edge
		if (fabs(P1[1] - P2[1]) < TOL){
			idLE = 1;
		}
		else if (fabs(P2[1] - P3[1]) < TOL){
			idLE = 2;
		}
		else if (fabs(P3[1] - P1[1]) < TOL){
			idLE = 3;
		}
		else{ // if no horizontal edge found, choose the longest edge at (pseudo-) random to break possible mesh symmetries and cycles
			idLE = rand() % 3 + 1; // random integer number 1, 2, 3
		}
	}
	else if (((fabs(e1 - e2) < TOL) && (e1 > e3)) ||
		((fabs(e1 - e3) < TOL) && (e1 > e2)) ||
		((fabs(e2 - e3) < TOL) && (e2 > e1))){ // two edges are of the same length

		// Take the longest edge at (pseudo-) random out of the two choices
		if (fabs(e1 - e2) < TOL){
			idLE = rand() % 2 + 1; // random integer number 1, 2
		}
		else if (fabs(e1 - e3) < TOL){
			idLE = rand() % 2 + 1; // random integer number 1, 3
			if (idLE == 2){
				idLE = 3;
			}
		}
		else if (fabs(e2 - e3) < TOL){
			idLE = rand() % 2 + 2; // random integer number 2, 3
		}
	}
	else{ // one unique longest edge exists
		if (e1 > e2){
			if (e1 > e3){
				idLE = 1;
			}
			else{
				idLE = 3;
			}
		}
		else{
			if (e2 > e3){
				idLE = 2;
			}
			else{
				idLE = 3;
			}
		}
	}

	// Test success
	if (idLE == -1){
		printf("LE in triangle %d not found.\n", id + 1);
	}

	return idLE;
}

// Find the LE-neighbour of id-th triangle
int longest_edge_neighbour(int** t, int* mark, int idLE, int Nt, int id){
	int idLENT = -1;
	int id1, id2, i1, i2, i3;

	// Extract appropriate edge indices
	if (idLE == 1){
		id1 = t[0][id];
		id2 = t[1][id];
	}
	else if (idLE == 2){
		id1 = t[1][id];
		id2 = t[2][id];
	}
	else if (idLE == 3){
		id1 = t[2][id];
		id2 = t[0][id];
	}

	// Find id1, id2 in prt - neighboring triangle has the same indices
	int found = 0;
	int i = 0;
	while (i < Nt && found == 0){
		if (i != id && mark[i] != -1){ // mark[i] == -1 means that i-th triangle is marked for deletion, skip it
			i1 = t[0][i];
			i2 = t[1][i];
			i3 = t[2][i];

			if (id1 == i1){
				if (id2 == i2 || id2 == i3){
					idLENT = i;
					found = 1;
				}
			}
			else if (id1 == i2){
				if (id2 == i1 || id2 == i3){
					idLENT = i;
					found = 1;
				}
			}
			else if (id1 == i3){
				if (id2 == i1 || id2 == i2){
					idLENT = i;
					found = 1;
				}
			}
		}

		// Proceed to the next triangle
		i++;
	}

	return idLENT;
}

// Find the LEPP of K-th triangle
void evaluate_LEPP(int** t, int Nt, int* mark, int K, int* LEPP){
	int next;
	int count = 0;
	int end = 0;

	LEPP[0] = 0;
	LEPP[count + 1] = K;
	while (end == 0){
		count++;
		next = longest_edge_neighbour(t, mark, 1, Nt, LEPP[count]); // MATLAB indexing
		if (next == -1){
			LEPP[count + 1] = 0; // LEPP ends on physical boundary
			end = 1;
		}
		else if (next == LEPP[count - 1] && count > 1){
			LEPP[count + 1] = 1; // LEPP ends by two triangles
			end = 1;
		}
		else{
			LEPP[count + 1] = next;
		}

		// Test LEPP buffer fullness
		if (count == NLEPP - 2){
			mexErrMsgTxt("Buffer is full, enlarge NLEPP and recompile.");
		}
	}
	LEPP[0] = count;
}
