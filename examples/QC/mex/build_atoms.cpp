// build_atoms: builds a database of all atoms.
// Implementation below exploits specific numbering of atoms from -SizeX to
// +SizeX and from -SizeY to +SizeY.

#include <stdio.h>
#include <math.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mex.h"
#include "matrix.h"

/******************** Structure declarations ******************/
struct ATOM{
	double R[2]; // coordinates of an atom in undeformed configuration, R = [x,y]
	int NeighbourList[9]; // length followed by a list of all the nearest neighbours, max length 9
};

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for one input argument
	if (nrhs != 1){
		mexErrMsgTxt("One input argument required.");
	}

	// Copy the input variables
	double *prrhs;
	prrhs = mxGetPr(prhs[0]);
	int Nrhs = (int)mxGetN(prhs[0]);

	// Test for the correct length of the input vector
	if (Nrhs != 4){
		mexErrMsgTxt("Input vector length must equal to 4.");
	}
	double SizeX = prrhs[0];
	double SizeY = prrhs[1];
	double dSizeX = prrhs[2];
	double dSizeY = prrhs[3];

	// Build database of all atoms
	int i, j;
	int NxAtoms = (int)floor(2 * SizeX / dSizeX) + 1; // number of atoms in x-direction
	int NyAtoms = (int)floor(2 * SizeY / dSizeY) + 1; // number of atoms in y-direction
	int NAtoms = NxAtoms*NyAtoms; // number of atoms
	ATOM *atoms = new ATOM[NAtoms]; // initialize the database of all atoms

	// Test for succesfull initialization
	if (atoms == NULL){
		mexErrMsgTxt("Not enough memory for atoms database.");
	}

	// Fill in the undeformed positions
	int count = 0;
	for (j = 0; j < NyAtoms; j++){
		for (i = 0; i < NxAtoms; i++){
			atoms[count].R[0] = -SizeX + i*dSizeX; // x-coordinate
			atoms[count].R[1] = -SizeY + j*dSizeY; // y-coordinate
			count++;
		}
	}

	// Build NeighbourList for all atoms
#pragma omp parallel for
	for (i = 0; i < NAtoms; i++){
		if (atoms[i].R[0] == -SizeX){
			if (atoms[i].R[1] == -SizeY){ // [-SizeX,-SizeY] corner
				atoms[i].NeighbourList[0] = 3; // 3 neighbours
				atoms[i].NeighbourList[1] = i + 2; // right, MATLAB indexing
				atoms[i].NeighbourList[2] = i + 1 + NxAtoms; // up
				atoms[i].NeighbourList[3] = i + 2 + NxAtoms; // up right
			}
			if (atoms[i].R[1] == SizeY){ // [-SizeX,SizeY] corner
				atoms[i].NeighbourList[0] = 3; // 3 neighbours
				atoms[i].NeighbourList[1] = i + 1 - NxAtoms; // down, MATLAB indexing
				atoms[i].NeighbourList[2] = i + 2 - NxAtoms; // down right
				atoms[i].NeighbourList[3] = i + 2; // right
			}
			if ((atoms[i].R[1] > -SizeY) && (atoms[i].R[1] < SizeY)){ // [-SizeX,Y] side
				atoms[i].NeighbourList[0] = 5; // 5 neighbours
				atoms[i].NeighbourList[1] = i + 1 - NxAtoms; // down, MATLAB indexing
				atoms[i].NeighbourList[2] = i + 2 - NxAtoms; // down right
				atoms[i].NeighbourList[3] = i + 2; // right
				atoms[i].NeighbourList[4] = i + 1 + NxAtoms; // up
				atoms[i].NeighbourList[5] = i + 2 + NxAtoms; // up right
			}
		}
		if (atoms[i].R[0] == SizeX){
			if (atoms[i].R[1] == -SizeY){ // [SizeX,-SizeY] corner
				atoms[i].NeighbourList[0] = 3; // 3 neighbours
				atoms[i].NeighbourList[1] = i; // left, MATLAB indexing
				atoms[i].NeighbourList[2] = i + NxAtoms; // up left
				atoms[i].NeighbourList[3] = i + 1 + NxAtoms; // up
			}
			if (atoms[i].R[1] == SizeY){ // [SizeX,SizeY] corner
				atoms[i].NeighbourList[0] = 3; // 3 neighbours
				atoms[i].NeighbourList[1] = i - NxAtoms; // down left, MATLAB indexing
				atoms[i].NeighbourList[2] = i + 1 - NxAtoms; // down
				atoms[i].NeighbourList[3] = i; // left
			}
			if ((atoms[i].R[1] > -SizeY) && (atoms[i].R[1] < SizeY)){ // [SizeX,Y] side
				atoms[i].NeighbourList[0] = 5; // 5 neighbours
				atoms[i].NeighbourList[1] = i - NxAtoms; // down left
				atoms[i].NeighbourList[2] = i + 1 - NxAtoms; // down
				atoms[i].NeighbourList[3] = i; // left, MATLAB indexing
				atoms[i].NeighbourList[4] = i + NxAtoms; // up left
				atoms[i].NeighbourList[5] = i + 1 + NxAtoms; // up
			}
		}
		if ((atoms[i].R[1] == -SizeY) && (atoms[i].R[0] > -SizeX) &&
			(atoms[i].R[0] < SizeX)){ // [X,-SizeY] side
			atoms[i].NeighbourList[0] = 5; // 5 neighbours
			atoms[i].NeighbourList[1] = i; // left, MATLAB indexing
			atoms[i].NeighbourList[2] = i + 2; // right
			atoms[i].NeighbourList[3] = i + NxAtoms; // up left
			atoms[i].NeighbourList[4] = i + 1 + NxAtoms; // up
			atoms[i].NeighbourList[5] = i + 2 + NxAtoms; // up right
		}
		if ((atoms[i].R[1] == SizeY) && (atoms[i].R[0] > -SizeX) &&
			(atoms[i].R[0] < SizeX)){ // [X,SizeY] side
			atoms[i].NeighbourList[0] = 5; // 5 neighbours
			atoms[i].NeighbourList[1] = i - NxAtoms; // down left, MATLAB indexing
			atoms[i].NeighbourList[2] = i + 1 - NxAtoms; // down
			atoms[i].NeighbourList[3] = i + 2 - NxAtoms; // down right
			atoms[i].NeighbourList[4] = i; // left
			atoms[i].NeighbourList[5] = i + 2; // right
		}
		if ((atoms[i].R[0] > -SizeX) && (atoms[i].R[0]<SizeX) &&
			(atoms[i].R[1]>-SizeY) && (atoms[i].R[1] < SizeY)){ // all internal atoms
			atoms[i].NeighbourList[0] = 8; // 8 neighbours
			atoms[i].NeighbourList[1] = i - NxAtoms; // down left
			atoms[i].NeighbourList[2] = i + 1 - NxAtoms; // down
			atoms[i].NeighbourList[3] = i + 2 - NxAtoms; // down right
			atoms[i].NeighbourList[4] = i; // left, MATLAB indexing
			atoms[i].NeighbourList[5] = i + 2; // right
			atoms[i].NeighbourList[6] = i + NxAtoms; // up left
			atoms[i].NeighbourList[7] = i + 1 + NxAtoms; // up
			atoms[i].NeighbourList[8] = i + 2 + NxAtoms; // up right
		}
	}

	// Send out the data back to MATLAB
	nlhs = 1;
	mwSize dims[] = { 1, NAtoms };
	const char *field_names[] = { "R", "NeighbourList", "BondList" };
	plhs[0] = mxCreateStructArray(2, dims, 3, field_names);
	double *prR, *prNeighbourList, *prBondList;
	mxArray *array_R, *array_NeighbourList, *array_BondList;
	int name_R = mxGetFieldNumber(plhs[0], "R");
	int name_NeighbourList = mxGetFieldNumber(plhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(plhs[0], "BondList");

	// Populate the output with computed data (mx functions are not thread-safe)
	for (i = 0; i < NAtoms; i++){

		// Create arrays
		array_R = mxCreateDoubleMatrix(1, 2, mxREAL);
		array_NeighbourList = mxCreateDoubleMatrix(1,
			atoms[i].NeighbourList[0], mxREAL);
		array_BondList = mxCreateDoubleMatrix(1,
			atoms[i].NeighbourList[0], mxREAL);

		// Get data arrays
		prR = mxGetPr(array_R);
		prNeighbourList = mxGetPr(array_NeighbourList);
		prBondList = mxGetPr(array_BondList);

		// Populate data arrays
		for (j = 0; j < 2; j++){
			prR[j] = atoms[i].R[j];
		}
		for (j = 0; j < atoms[i].NeighbourList[0]; j++){
			prNeighbourList[j] = atoms[i].NeighbourList[j + 1];
		}

		// Assign arrays to output structure
		mxSetFieldByNumber(plhs[0], i, name_R, array_R);
		mxSetFieldByNumber(plhs[0], i, name_NeighbourList,
			array_NeighbourList);
		mxSetFieldByNumber(plhs[0], i, name_BondList, array_BondList);
	}

	// Delete atoms database
	delete[] atoms;
}
