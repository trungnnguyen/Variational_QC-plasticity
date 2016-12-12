// build_bonds: builds a database of all bonds, updates BondLists of atoms'
// database.
// Implementation below exploits specific numbering of atoms from -SizeX to
// +SizeX and from -SizeY to +SizeY.

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mex.h"
#include "matrix.h"

/******************** Structure declarations ******************/
struct BOND{
	int ID; // ID of a bond
	int Atoms[2]; // atoms connected by this bond
};

/******************** Function declarations *******************/
int distID(const void *BondA, const void *BondB); // comparison function w.r.t. ID
int distAtom1(const void *BondA, const void *BondB); // comparison function w.r.t. the first atom
int distAtom2(const void *BondA, const void *BondB); // comparison function w.r.t. the second atom

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for six input arguments
	if (nrhs != 6){
		mexErrMsgTxt("Six input arguments required.");
	}

	// Copy input variables
	double* prRect = mxGetPr(prhs[1]);
	if ((int)mxGetN(prhs[1]) != 6){
		mexErrMsgTxt("Input vector specifying a rectangle must be of length 6.");
	}
	double SizeX = prRect[0];
	double SizeY = prRect[1];
	double dSizeX = prRect[2];
	double dSizeY = prRect[3];
	double RigidX = prRect[4];
	double RigidY = prRect[5];
	int NAtoms = (int)mxGetNumberOfElements(prhs[0]);
	double* PotentialI = mxGetPr(prhs[2]);
	double* PotentialM = mxGetPr(prhs[3]);
	if ((int)mxGetN(prhs[2]) != 4 || (int)mxGetN(prhs[3]) != 4){
		mexErrMsgTxt("Input vector specifying potentials must be of length 4.");
	}
	double *prSW = mxGetPr(prhs[4]);
	int SW = (int)prSW[0]; // switch for different locations of the inclusion
	double* prTOL = mxGetPr(prhs[5]);
	double TOL = prTOL[0]; // a positive distance < TOL is treated as zero

	// Initialize the data
	int i, j;
	int NxAtoms = (int)floor(2 * SizeX / dSizeX) + 1; // number of atoms in x-direction
	int NyAtoms = (int)floor(2 * SizeY / dSizeY) + 1; // number of atoms in y-direction
	int NBonds = NyAtoms*(NxAtoms - 1) + NxAtoms*(NyAtoms - 1) +
		2 * (NyAtoms - 1)*(NxAtoms - 1);
	BOND* BondsA = new BOND[NBonds]; // bonds database
	BOND* BondsB = new BOND[NBonds]; // working bonds database
	if (BondsA == NULL || BondsB == NULL){
		mexErrMsgTxt("Not enough memory for Bonds.");
	}

	// Build BondsA and BondsB databases
	// Horizontal bonds
	int counter = 0;
	for (j = 0; j < NyAtoms; j++){
		for (i = 0; i < (NxAtoms - 1); i++){
			BondsA[counter].ID = counter + 1; // MATLAB indexing
			BondsA[counter].Atoms[0] = i + 1 + j*NxAtoms;
			BondsA[counter].Atoms[1] = i + 2 + j*NxAtoms;
			counter++;
		}
	}

	// Vertical bonds
	for (j = 0; j < NxAtoms; j++){
		for (i = 0; i < (NyAtoms - 1); i++){
			BondsA[counter].ID = counter + 1; // MATLAB indexing
			BondsA[counter].Atoms[0] = j + 1 + i*NxAtoms;
			BondsA[counter].Atoms[1] = j + 1 + (i + 1)*NxAtoms;
			counter++;
		}
	}

	// Ascending diagonals
	for (j = 0; j < (NyAtoms - 1); j++){
		for (i = 0; i < (NxAtoms - 1); i++){
			BondsA[counter].ID = counter + 1; // MATLAB indexing
			BondsA[counter].Atoms[0] = i + 1 + j*NxAtoms;
			BondsA[counter].Atoms[1] = i + NxAtoms + 2 + j*NxAtoms;
			counter++;
		}
	}

	// Descending diagonals
	for (j = 0; j < (NyAtoms - 1); j++){
		for (i = 0; i < (NxAtoms - 1); i++){
			BondsA[counter].ID = counter + 1; // MATLAB indexing
			BondsA[counter].Atoms[0] = i + 1 + NxAtoms + j*NxAtoms;
			BondsA[counter].Atoms[1] = i + 2 + j*NxAtoms;
			counter++;
		}
	}

	// Copy BondsA to BondsB
	for (i = 0; i < NBonds; i++){
		BondsB[i].ID = BondsA[i].ID;
		BondsB[i].Atoms[0] = BondsA[i].Atoms[0];
		BondsB[i].Atoms[1] = BondsA[i].Atoms[1];
	}

	// Sort both databases w.r.t. first and secod atoms
#pragma omp parallel sections
	{
#pragma omp section
		{qsort(BondsA, NBonds, sizeof(BOND), distAtom1); } // sort w.r.t. the first atom ID
#pragma omp section
		{qsort(BondsB, NBonds, sizeof(BOND), distAtom2); } // sort w.r.t. the second atom ID
	}

	// Create vectors of multiplicity
	int* multA = new int[NAtoms];
	int m = 0;
	for (i = 0; i < NAtoms; i++){ // for first atoms in BondsA
		j = 0;
		while (BondsA[m].Atoms[0] == (i + 1)){ // MATLAB indexing
			j++;
			m++;
		}
		multA[i] = j;
	}
	int* multB = new int[NAtoms];
	m = 0;
	for (i = 0; i < NAtoms; i++){ // for first atoms in BondsB
		j = 0;
		while (BondsB[m].Atoms[1] == (i + 1)){ // MATLAB indexing
			j++;
			m++;
		}
		multB[i] = j;
	}

	// Create vectors of cummulative sums
	int* cumsumA = new int[NAtoms];
	int* cumsumB = new int[NAtoms];
	cumsumA[0] = multA[0];
	cumsumB[0] = multB[0];
	for (i = 1; i < NAtoms; i++){
		cumsumA[i] = cumsumA[i - 1] + multA[i];
		cumsumB[i] = cumsumB[i - 1] + multB[i];
	}

	// For all atoms and all neighbours, find ID of their bonds and store the data
	double *NeighbourList, *BondList;
	bool foundA, foundB;
	int alpha, beta, tempMin, tempMax;
	mxArray *pBondList;
	int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	for (i = 0; i < NAtoms; i++){ // for cycle over all atoms
		alpha = i + 1; // ID of an atom alpha

		// Get data for atom alpha
		NeighbourList = mxGetPr(mxGetFieldByNumber(prhs[0], alpha - 1,
			name_NeighbourList));
		pBondList = mxGetFieldByNumber(prhs[0], alpha - 1, name_BondList);
		BondList = mxGetPr(pBondList);
		for (j = 0; j < mxGetN(pBondList); j++){ // for cycle over all nearest neighbours
			beta = (int)NeighbourList[j]; // ID of atom beta, MATLAB indexing

			// Find in BondsA, BondsB ID of the bond among atoms alpha and beta;
			// we know that two atoms are connected only by one bond
			foundA = false;
			foundB = false;
			counter = 0;
			if (alpha == 1){
				tempMin = 0;
				tempMax = cumsumA[(alpha - 1)];
			}
			else{
				tempMin = cumsumA[(alpha - 1) - 1];
				tempMax = cumsumA[(alpha - 1)];
			}
			for (m = tempMin; m < tempMax; m++){ // MATLAB indexing
				if (BondsA[m].Atoms[1] == beta){
					foundA = true;
					counter = m;
					break;
				}
			}
			if (foundA == false){ // if not found, then reverse the order beta and alpha
				if (alpha == 1){
					tempMin = 0;
					tempMax = cumsumB[(alpha - 1)];
				}
				else{
					tempMin = cumsumB[(alpha - 1) - 1];
					tempMax = cumsumB[(alpha - 1)];
				}
				for (m = tempMin; m < tempMax; m++){ // MATLAB indexing
					if (BondsB[m].Atoms[0] == beta){
						foundB = true;
						counter = m;
						break;
					}
				}
			}
			if (foundA == true){
				BondList[j] = (double)BondsA[counter].ID; // store found ID to BondList
			}
			else if (foundB == true){
				BondList[j] = (double)BondsB[counter].ID; // store found ID to BondList
			}
			else{
				printf("Unable to find bond for atoms %u and %u.\n",
					alpha, beta);
			}
		}
	}
	delete[] multA;
	delete[] multB;
	delete[] cumsumA;
	delete[] cumsumB;
	delete[] BondsB;

	// Sort resulting BondsA database by IDs
	qsort(BondsA, NBonds, sizeof(BOND), distID);

	// Construct bonds database and send out the data back to MATLAB
	nlhs = 1;
	mwSize dims[] = { 1, NBonds };
	const char *field_names[] = { "Atoms", "Potential" };
	plhs[0] = mxCreateStructArray(2, dims, 2, field_names);

	// Populate output with computed data (mx functions are not thread-safe)
	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_Atoms = mxGetFieldNumber(plhs[0], "Atoms");
	int name_Potential = mxGetFieldNumber(plhs[0], "Potential");
	mxArray *array_Atoms, *array_Potential;
	double *prAtoms, *prPotential, *Ra, *Rb;
	for (i = 0; i < NBonds; i++){

		// Create arrays
		array_Atoms = mxCreateDoubleMatrix(1, 2, mxREAL);
		array_Potential = mxCreateDoubleMatrix(1, 4, mxREAL);

		// Get data arrays
		prAtoms = mxGetPr(array_Atoms);
		prPotential = mxGetPr(array_Potential);

		// And populate them
		for (j = 0; j < 2; j++){
			prAtoms[j] = BondsA[i].Atoms[j];
		}

		/*<<<<<<<<<<<<<<<<<<<<<<<< Insert rigid region >>>>>>>>>>>>>>>>>>>>>>>>*/
		// Get data for atoms alpha and beta
		alpha = BondsA[i].Atoms[0]; // bond's atom alpha, MATLAB indexing
		beta = BondsA[i].Atoms[1]; // bond's atom beta, MATLAB indexing
		Ra = mxGetPr(mxGetFieldByNumber(prhs[0], alpha - 1,
			name_R)); // initial coordinates of an atom alpha
		Rb = mxGetPr(mxGetFieldByNumber(prhs[0], beta - 1,
			name_R)); // initial coordinates of an atom beta

		switch (SW){

			/*<<<<<<<<<<<<<<<<<<<<<<<< Uniform loading case >>>>>>>>>>>>>>>>>>>>>>>>*/
		case 0:

			// If atoms alpha and beta are within the central subdomain
			if ((fabs(Ra[0]) < RigidX + TOL) && (fabs(Ra[1]) < RigidY + TOL)
				&& (fabs(Rb[0]) < RigidX + TOL) && (fabs(Rb[1]) < RigidY + TOL)){

				// Assign parameters stored in PotentialI
				for (j = 0; j < 4; j++){
					prPotential[j] = PotentialI[j];
				}
			}

			// If alpha and beta are not within the central subdomain
			else{

				// Assign parameters stored in PotentialM
				for (j = 0; j < 4; j++){
					prPotential[j] = PotentialM[j];
				}
			}

			// Assign half Young's modulus to bonds on the physical boundary
			// If atoms alpha and beta lie on the physical boundary
			if (((fabs(fabs(Ra[0]) - SizeX) < TOL) &&
				(fabs(fabs(Rb[0]) - SizeX) < TOL)) ||
				((fabs(fabs(Ra[1]) - SizeY) < TOL) &&
				(fabs(fabs(Rb[1]) - SizeY) < TOL))){

				// Divide the Young modulus and initial yield stress by two: corresponds to halving the sectional area A
				prPotential[0] = 0.5*PotentialM[0]; // halve Young's modulus
				prPotential[2] = 0.5*PotentialM[2]; // halve yield stress \sigma_0
			}
			break;

			/*<<<<<<<<<<<<<<<<<<<<<<<< Pure bending case >>>>>>>>>>>>>>>>>>>>>>>>*/
		case 1:

			// If atoms alpha and beta are within the bottom-edge subdomain
			if ((fabs(Ra[0]) < RigidX + TOL) &&
				(Ra[1] < RigidY - SizeY + TOL) &&
				(fabs(Rb[0]) < RigidX + TOL) &&
				(Rb[1] < RigidY - SizeY + TOL)){

				// Assign parameters stored in PotentialI
				for (j = 0; j < 4; j++){
					prPotential[j] = PotentialI[j];
				}
			}

			// If alpha and beta are not within the bottom-edge subdomain
			else{

				// Assign parameters stored in PotentialM
				for (j = 0; j < 4; j++){
					prPotential[j] = PotentialM[j];
				}
			}

			// Assign half Young's modulus to bonds on the \Gamma_2 and \Gamma_4 parts of the boundary
			// If atoms alpha and beta lie on the \Gamma_2 or \Gamma_4 parts of the boundary
			if ((fabs(fabs(Ra[0]) - SizeX) < TOL) &&
				(fabs(fabs(Rb[0]) - SizeX) < TOL)){

				// Divide the Young modulus and initial yield stress by two: corresponds to halving the sectional area A
				prPotential[0] = 0.5*PotentialM[0]; // halve Young's modulus
				prPotential[2] = 0.5*PotentialM[2]; // halve yield stress \sigma_0
			}
			break;

		case 2:

			// Assign parameters stored in PotentialM
			for (j = 0; j < 4; j++){
				prPotential[j] = PotentialM[j];
			}
			break;
		}

		// Assign output arrays
		mxSetFieldByNumber(plhs[0], i, name_Atoms, array_Atoms);
		mxSetFieldByNumber(plhs[0], i, name_Potential, array_Potential);
	}
	delete[] BondsA;
}

/******************** Function definitions ********************/
// Comparison function w.r.t. ID
int distID(const void *BondA, const void *BondB){
	BOND *A = (BOND*)BondA,
		*B = (BOND*)BondB;
	double r;
	r = A->ID - B->ID;
	if (r < 0)
		return -1;
	else if (r > 0)
		return 1;
	else
		return 0;
}

// Comparison function w.r.t. first atom
int distAtom1(const void *BondA, const void *BondB){
	BOND *A = (BOND*)BondA,
		*B = (BOND*)BondB;
	double r;
	r = A->Atoms[0] - B->Atoms[0];
	if (r < 0)
		return -1;
	else if (r > 0)
		return 1;
	else
		return 0;
}

// Comparison function w.r.t. second atom
int distAtom2(const void *BondA, const void *BondB){
	BOND *A = (BOND*)BondA,
		*B = (BOND*)BondB;
	double r;
	r = A->Atoms[1] - B->Atoms[1];
	if (r < 0)
		return -1;
	else if (r > 0)
		return 1;
	else
		return 0;
}
