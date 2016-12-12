// sort_atoms_QC: builds a database of triangles, vector of repatoms,
// and I, J, S arrays for the COO representation of the matrix Phi.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mex.h"
#include "matrix.h"

/******************* Definition of constants ******************/
const int FLAG = 2; // definition of triangle's center: 1 - centroid, 2 - incenter, 3 - circumcenter

/******************* Declaration of global variables ******************/
double TOL; // a positive number < TOL is considered as zero

/******************** Structure declarations ******************/
struct TRIANGLE{
	double P1[2], P2[2], P3[2]; // triangle vertex coordinates
	double T[2]; // center of a triangle: centroid, circumcenter
	double r2; // squared radius of circumcircle to a triangle
	double Area; // area of a triangle
	int *IntAtoms, *EdgeAtoms, *VertexAtoms; // IDs of appropriate atoms
	double *Lint, *Ledge; // vectors storing barycentric coordinates for appropriate atoms
	int *NeighTriangles; // neighbouring triangles followed by their IDs
};

struct ATOM{
	double P[2]; // coordinates of an atom
	double Dist2; // square of the distance from a point
};

/******************** Function declarations *******************/
void TriagCentroid(double Cg[], double *P1, double *P2, double *P3, int flag); // according to flag computes internal point of a triangle
int distInt(const void *a, const void *b); // comparison function for integers

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for four input arguments
	if (nrhs != 4){
		mexErrMsgTxt("Four input arguments required.");
	}

	// Copy input variables
	double* prp = mxGetPr(prhs[0]);
	int Mp = (int)mxGetM(prhs[0]);
	int Np = (int)mxGetN(prhs[0]);
	double* prt = mxGetPr(prhs[1]);
	int Mt = (int)mxGetM(prhs[1]);
	int Nt = (int)mxGetN(prhs[1]);
	int NAtoms = (int)mxGetNumberOfElements(prhs[2]);
	double* prTOL = mxGetPr(prhs[3]);
	TOL = prTOL[0]; // a positive number < TOL is treated as zero

	// Allocate memory
	TRIANGLE* Triangles = new TRIANGLE[Nt];
	ATOM* Atoms = new ATOM[NAtoms];
	if (Triangles == NULL || Atoms == NULL){
		mexErrMsgTxt("Not enough memory for databases.");
	}

	// Fill up Atoms
	int name_R = mxGetFieldNumber(prhs[2], "R");
	double* Ra;
	int i;
	for (i = 0; i < NAtoms; i++){
		Ra = mxGetPr(mxGetFieldByNumber(prhs[2], i, name_R));
		Atoms[i].P[0] = Ra[0]; // x coordinate
		Atoms[i].P[1] = Ra[1]; // y coordinate
	}

	// Fill up Triangles
	double r[3], Cg[2];
	for (i = 0; i < Nt; i++){
		Triangles[i].P1[0] = prp[(int)(prt[Mt*i] - 1)*Mp]; // x coordinate of P1
		Triangles[i].P1[1] = prp[(int)(prt[Mt*i] - 1)*Mp + 1]; // y coordinate of P1
		Triangles[i].P2[0] = prp[(int)(prt[1 + Mt*i] - 1)*Mp];
		Triangles[i].P2[1] = prp[(int)(prt[1 + Mt*i] - 1)*Mp + 1];
		Triangles[i].P3[0] = prp[(int)(prt[2 + Mt*i] - 1)*Mp];
		Triangles[i].P3[1] = prp[(int)(prt[2 + Mt*i] - 1)*Mp + 1];
		TriagCentroid(Cg, Triangles[i].P1, Triangles[i].P2,
			Triangles[i].P3, FLAG); // FLAG: 1 - centroid, 2 - incenter, 3 - circumcenter
		Triangles[i].T[0] = Cg[0];
		Triangles[i].T[1] = Cg[1];
		Triangles[i].Area = fabs((Triangles[i].P1[0] *
			(Triangles[i].P2[1] - Triangles[i].P3[1]) +
			Triangles[i].P2[0] * (Triangles[i].P3[1] - Triangles[i].P1[1]) +
			Triangles[i].P3[0] * (Triangles[i].P1[1] - Triangles[i].P2[1])) / 2);

		// Squared radius of circumcircle with center in T
		r[0] = pow(Cg[0] - Triangles[i].P1[0], 2) +
			pow(Cg[1] - Triangles[i].P1[1], 2);
		r[1] = pow(Cg[0] - Triangles[i].P2[0], 2) +
			pow(Cg[1] - Triangles[i].P2[1], 2);
		r[2] = pow(Cg[0] - Triangles[i].P3[0], 2) +
			pow(Cg[1] - Triangles[i].P3[1], 2);
		Triangles[i].r2 = std::max(std::max(r[0], r[1]), r[2]);
	}





	/*<<<<<<<<<<<<<<<<<<<<<<<< Sort atoms within triangles >>>>>>>>>>>>>>>>>>>>>>>>*/
	// Find atoms inside triangles
	int Phicount = 0; // number of nnz in Phi matrix including duplicities
#pragma omp parallel
	{
		double _alpha, _beta, _gamma, *_Ra; // barycentric coordinates and atom's position vector
		int _Phicount = 0; // number of nnz in Phi matrix including duplicities
		int* _IntAtoms = new int[NAtoms]; // temporary fields for storing atoms within a triangle i
		int* _EdgeAtoms = new int[NAtoms];
		int* _VertexAtoms = new int[NAtoms];
		double *_Lint = new double[3 * NAtoms]; // temporary fields for storing barycentric coordinates for atoms within a triangle i
		double *_Ledge = new double[3 * NAtoms];

		// Loop over all triangles
#pragma omp for 
		for (i = 0; i < Nt; i++){
			int _nInt = 0;
			int _nEdge = 0;
			int _nVertex = 0;

			// Find which atoms are inside i-th triangle through barycentric coordinates
			for (int _j = 0; _j < NAtoms; _j++){

				if ((pow(Triangles[i].T[0] - Atoms[_j].P[0], 2)
					+ pow(Triangles[i].T[1] - Atoms[_j].P[1], 2)) <
					Triangles[i].r2 + TOL){ // if atom is inside circumcircle

					// Compute barycentric coordinates
					_alpha = ((Triangles[i].P2[1] - Triangles[i].P3[1])*
						(Atoms[_j].P[0] - Triangles[i].P3[0]) +
						(Triangles[i].P3[0] - Triangles[i].P2[0])*
						(Atoms[_j].P[1] - Triangles[i].P3[1])) /
						((Triangles[i].P2[1] - Triangles[i].P3[1])*
						(Triangles[i].P1[0] - Triangles[i].P3[0]) +
						(Triangles[i].P3[0] - Triangles[i].P2[0])*
						(Triangles[i].P1[1] - Triangles[i].P3[1]));
					_beta = ((Triangles[i].P3[1] - Triangles[i].P1[1])*
						(Atoms[_j].P[0] - Triangles[i].P3[0]) +
						(Triangles[i].P1[0] - Triangles[i].P3[0])*
						(Atoms[_j].P[1] - Triangles[i].P3[1])) /
						((Triangles[i].P2[1] - Triangles[i].P3[1])*
						(Triangles[i].P1[0] - Triangles[i].P3[0]) +
						(Triangles[i].P3[0] - Triangles[i].P2[0])*
						(Triangles[i].P1[1] - Triangles[i].P3[1]));
					_gamma = 1 - _alpha - _beta;

					// Coordinates have to be in [0,1] in order to be inside the triangle
					if ((_alpha >= -TOL) && (_beta >= -TOL) && (_gamma >= -TOL)){

						// Decide where the atom lies: inside, edge, vertex
						if ((fabs(_alpha) > TOL) && (fabs(_beta) > TOL) &&
							(fabs(_gamma) > TOL)){
							_IntAtoms[_nInt] = _j + 1; // MATLAB indexing
							_Lint[3 * _nInt + 0] = _alpha;
							_Lint[3 * _nInt + 1] = _beta;
							_Lint[3 * _nInt + 2] = _gamma;
							_nInt++;
						}
						else if (((fabs(_alpha) < TOL) && (fabs(_beta) > TOL) &&
							(fabs(_gamma) > TOL)) || ((fabs(_alpha) > TOL) &&
							(fabs(_beta) < TOL) && (fabs(_gamma) > TOL)) ||
							((fabs(_alpha) > TOL) && (fabs(_beta) > TOL) &&
							(fabs(_gamma) < TOL))){
							_EdgeAtoms[_nEdge] = _j + 1; // MATLAB indexing
							_Ledge[3 * _nEdge + 0] = _alpha;
							_Ledge[3 * _nEdge + 1] = _beta;
							_Ledge[3 * _nEdge + 2] = _gamma;
							_nEdge++;
						}
						else if (((fabs(_alpha) < TOL) && (fabs(_beta) < TOL) &&
							(fabs(_gamma) > TOL)) || ((fabs(_alpha) < TOL) &&
							(fabs(_beta) > TOL) && (fabs(_gamma) < TOL)) ||
							((fabs(_alpha) > TOL) && (fabs(_beta) < TOL) &&
							(fabs(_gamma) < TOL))){
							_VertexAtoms[_nVertex] = _j + 1; // MATLAB indexing
							_nVertex++;
						}
						else{
							printf("\nAtom %u not assigned properly in triangle %u\n.", _j, i + 1);
						}
					}
				}
			}

			// Update the list of atoms for i-th triangle
			Triangles[i].IntAtoms = new int[_nInt + 1];
			Triangles[i].EdgeAtoms = new int[_nEdge + 1];
			Triangles[i].VertexAtoms = new int[_nVertex + 1];
			Triangles[i].IntAtoms[0] = _nInt;
			Triangles[i].EdgeAtoms[0] = _nEdge;
			Triangles[i].VertexAtoms[0] = _nVertex;
			Triangles[i].Lint = new double[3 * _nInt];
			Triangles[i].Ledge = new double[3 * _nEdge];
			for (int _k = 0; _k < _nInt; _k++){
				Triangles[i].IntAtoms[_k + 1] = _IntAtoms[_k];
			}
			for (int _k = 0; _k < _nEdge; _k++){
				Triangles[i].EdgeAtoms[_k + 1] = _EdgeAtoms[_k];
			}

			// Store VertexAtoms according to P1, P2, and P3 coordinates rather than by distance from T
			if (_nVertex == 3){ // if three vertices are found, sort them according to P1, P2, and P3
				Triangles[i].VertexAtoms[0] = 3;
				for (int _k = 0; _k < 3; _k++){
					_Ra = mxGetPr(mxGetFieldByNumber(prhs[2],
						_VertexAtoms[_k] - 1, name_R));
					if ((fabs(Triangles[i].P1[0] - _Ra[0]) < TOL) &&
						(fabs(Triangles[i].P1[1] - _Ra[1]) < TOL)){
						Triangles[i].VertexAtoms[1] = _VertexAtoms[_k];
					}
					else if ((fabs(Triangles[i].P2[0] - _Ra[0]) < TOL) &&
						(fabs(Triangles[i].P2[1] - _Ra[1]) < TOL)){
						Triangles[i].VertexAtoms[2] = _VertexAtoms[_k];
					}
					else if ((fabs(Triangles[i].P3[0] - _Ra[0]) < TOL) &&
						(fabs(Triangles[i].P3[1] - _Ra[1]) < TOL)){
						Triangles[i].VertexAtoms[3] = _VertexAtoms[_k];
					}
					else {
						printf("Vertex atom %d not allocated.\n", _VertexAtoms[_k]);
					}
				}
			}
			else {
				printf("Number of vertex atoms in triangle %u does not sum to 3.\n", i + 1);
				for (int _k = 0; _k < _nVertex; _k++){ // if three vertices not found, order does not matter
					Triangles[i].VertexAtoms[_k + 1] = _VertexAtoms[_k];
				}
			}

			// Store barycentric coordinates
			for (int _k = 0; _k < _nInt; _k++){
				Triangles[i].Lint[3 * _k + 0] = _Lint[3 * _k + 0];
				Triangles[i].Lint[3 * _k + 1] = _Lint[3 * _k + 1];
				Triangles[i].Lint[3 * _k + 2] = _Lint[3 * _k + 2];
			}
			for (int _k = 0; _k < _nEdge; _k++){
				Triangles[i].Ledge[3 * _k + 0] = _Ledge[3 * _k + 0];
				Triangles[i].Ledge[3 * _k + 1] = _Ledge[3 * _k + 1];
				Triangles[i].Ledge[3 * _k + 2] = _Ledge[3 * _k + 2];
			}

			// Update Phicount
#pragma omp atomic
			Phicount += 3 * (_nInt + _nEdge) + _nVertex;
		}

		// Delete temporary arrays
		delete[] _IntAtoms;
		delete[] _EdgeAtoms;
		delete[] _VertexAtoms;
		delete[] _Lint;
		delete[] _Ledge;
	}

	// Delete Atoms database
	delete[] Atoms;






	/*<<<<<<<<<<<<<<<<<<<<<<<< Construct triangle neighbours and repatoms >>>>>>>>>>>>>>>>>>>>>>>>*/
	// For each triangle construct a list of neighbouring triangles
#pragma omp parallel
	{
		bool _found;
		int _j, _TrNeigh, _Va, _Vb, _Vc, _Vd, _Ve, _Vf, _tempList[100];

#pragma omp for
		for (i = 0; i < Nt; i++){ // loop over all triangles
			_TrNeigh = 0; // number of neighbours found
			if (Triangles[i].VertexAtoms[0] == 3){
				_Va = Triangles[i].VertexAtoms[1]; // vertex atoms of i-th triangle, MATLAB indexing
				_Vb = Triangles[i].VertexAtoms[2];
				_Vc = Triangles[i].VertexAtoms[3];
				for (_j = 0; _j < Nt; _j++){
					if (_j != i){
						_found = 0;
						_Vd = Triangles[_j].VertexAtoms[1]; // vertex atoms of j-th triangle, MATLAB indexing
						_Ve = Triangles[_j].VertexAtoms[2];
						_Vf = Triangles[_j].VertexAtoms[3];
						if (((_Va == _Vd) && (_Vb == _Ve || _Vb == _Vf)) ||
							((_Va == _Ve) && (_Vb == _Vf)) ||
							((_Vb == _Vd) && (_Va == _Ve || _Va == _Vf)) ||
							((_Vb == _Ve) && (_Va == _Vf))){
							_found = 1;
							_tempList[_TrNeigh] = _j + 1; // MATLAB indexing
							_TrNeigh++;
						}
						if (_found == 0){
							if (((_Va == _Vd) && (_Vc == _Ve || _Vc == _Vf)) ||
								((_Va == _Ve) && (_Vc == _Vf)) ||
								((_Vc == _Vd) && (_Va == _Ve || _Va == _Vf)) ||
								((_Vc == _Ve) && (_Va == _Vf))){
								_found = 1;
								_tempList[_TrNeigh] = _j + 1; // MATLAB indexing
								_TrNeigh++;
							}
						}
						if (_found == 0){
							if (((_Vb == _Vd) && (_Vc == _Ve || _Vc == _Vf)) ||
								((_Vb == _Ve) && (_Vc == _Vf)) ||
								((_Vc == _Vd) && (_Vb == _Ve || _Vb == _Vf)) ||
								((_Vc == _Ve) && (_Vb == _Vf))){
								_found = 1;
								_tempList[_TrNeigh] = _j + 1; // MATLAB indexing
								_TrNeigh++;
							}
						}
					}
				}

				// Update NeighTriangles list of i-th triangle
				Triangles[i].NeighTriangles = new int[_TrNeigh + 1];
				Triangles[i].NeighTriangles[0] = _TrNeigh;
				for (_j = 0; _j < _TrNeigh; _j++){
					Triangles[i].NeighTriangles[_j + 1] = _tempList[_j];
				}
				if (_TrNeigh == 0){
					printf("No neighbour of triangle %u found.\n", i + 1);
				}
				if (_TrNeigh > 3){
					printf("%u neighbours of triangle %u found.\n", _TrNeigh, i + 1);
				}
			}
		}
	}

	// Construct triangles database and send out the data back to MATLAB
	nlhs = 5;
	mwSize dims[] = { 1, Nt };
	const char *field_names[] = { "P1", "P2", "P3", "T", "IntAtoms",
		"EdgeAtoms", "VertexAtoms", "NeighTriangles" };
	plhs[0] = mxCreateStructArray(2, dims, 8, field_names);

	// Populate output with computed data
	int j, m;
	int name_P1 = mxGetFieldNumber(plhs[0], "P1");
	int name_P2 = mxGetFieldNumber(plhs[0], "P2");
	int name_P3 = mxGetFieldNumber(plhs[0], "P3");
	int name_T = mxGetFieldNumber(plhs[0], "T");
	int name_IntAtoms = mxGetFieldNumber(plhs[0], "IntAtoms");
	int name_EdgeAtoms = mxGetFieldNumber(plhs[0], "EdgeAtoms");
	int name_VertexAtoms = mxGetFieldNumber(plhs[0], "VertexAtoms");
	int name_NeighTriangles = mxGetFieldNumber(plhs[0], "NeighTriangles");
	for (i = 0; i < Nt; i++){

		// Create arrays
		mxArray* array_P1 = mxCreateDoubleMatrix(1, 2, mxREAL);
		mxArray* array_P2 = mxCreateDoubleMatrix(1, 2, mxREAL);
		mxArray* array_P3 = mxCreateDoubleMatrix(1, 2, mxREAL);
		mxArray* array_T = mxCreateDoubleMatrix(1, 2, mxREAL);
		mxArray* array_IntAtoms = mxCreateDoubleMatrix(1,
			Triangles[i].IntAtoms[0], mxREAL);
		mxArray* array_EdgeAtoms = mxCreateDoubleMatrix(1,
			Triangles[i].EdgeAtoms[0], mxREAL);
		mxArray* array_VertexAtoms = mxCreateDoubleMatrix(1,
			Triangles[i].VertexAtoms[0], mxREAL);
		mxArray* array_NeighTriangles = mxCreateDoubleMatrix(1,
			Triangles[i].NeighTriangles[0], mxREAL);

		// Get data arrays
		double* prP1 = mxGetPr(array_P1);
		double* prP2 = mxGetPr(array_P2);
		double* prP3 = mxGetPr(array_P3);
		double* prT = mxGetPr(array_T);
		double* prIntAtoms = mxGetPr(array_IntAtoms);
		double* prEdgeAtoms = mxGetPr(array_EdgeAtoms);
		double* prVertexAtoms = mxGetPr(array_VertexAtoms);
		double* prNeighTriangles = mxGetPr(array_NeighTriangles);

		// Populate data arrays
		for (m = 0; m < 2; m++){
			prP1[m] = Triangles[i].P1[m];
		}
		for (m = 0; m < 2; m++){
			prP2[m] = Triangles[i].P2[m];
		}
		for (m = 0; m < 2; m++){
			prP3[m] = Triangles[i].P3[m];
		}
		for (m = 0; m < 2; m++){
			prT[m] = Triangles[i].T[m];
		}
		for (m = 0; m < Triangles[i].IntAtoms[0]; m++){
			prIntAtoms[m] = (double)Triangles[i].IntAtoms[m + 1];
		}
		for (m = 0; m < Triangles[i].EdgeAtoms[0]; m++){
			prEdgeAtoms[m] = (double)Triangles[i].EdgeAtoms[m + 1];
		}
		for (m = 0; m < Triangles[i].VertexAtoms[0]; m++){
			prVertexAtoms[m] = (double)Triangles[i].VertexAtoms[m + 1];
		}
		for (m = 0; m < Triangles[i].NeighTriangles[0]; m++){
			prNeighTriangles[m] = (double)Triangles[i].NeighTriangles[m + 1];
		}

		// Assign data arrays
		mxSetFieldByNumber(plhs[0], i, name_P1, array_P1);
		mxSetFieldByNumber(plhs[0], i, name_P2, array_P2);
		mxSetFieldByNumber(plhs[0], i, name_P3, array_P3);
		mxSetFieldByNumber(plhs[0], i, name_T, array_T);
		mxSetFieldByNumber(plhs[0], i, name_IntAtoms, array_IntAtoms);
		mxSetFieldByNumber(plhs[0], i, name_EdgeAtoms, array_EdgeAtoms);
		mxSetFieldByNumber(plhs[0], i, name_VertexAtoms, array_VertexAtoms);
		mxSetFieldByNumber(plhs[0], i, name_NeighTriangles, array_NeighTriangles);
	}

	// Construct Nodes array
	int* Nodes = new int[Nt * 3];
	for (i = 0; i < Nt; i++){
		for (m = 0; m < Triangles[i].VertexAtoms[0]; m++){
			Nodes[3 * i + m] = Triangles[i].VertexAtoms[m + 1];
		}
	}
	qsort(Nodes, 3 * Nt, sizeof(int), distInt);

	// Test for multiple ID occurences and zeros in Nodes vector
	int current, previous;
	i = 0;
	if (Nodes[0] == 0){
		while (Nodes[i] == 0 && i < 3 * Nt){
			i++;
		}
		printf("Possible loss of data in 2nd output, %u zeros reported.\n", i);
	}
	previous = Nodes[i];
	int counter = 0;
	for (j = i; j < 3 * Nt; j++){
		current = Nodes[j];
		if (current != previous){
			counter++;
			previous = current;
		}
	}

	// Populate AllNodes output, i.e. repatoms
	double *prAllNodes;
	int Nrepatoms = counter + 1;
	plhs[1] = mxCreateDoubleMatrix(Nrepatoms, 1, mxREAL);
	prAllNodes = mxGetPr(plhs[1]);
	previous = Nodes[i];
	prAllNodes[0] = previous;
	counter = 0;
	for (j = i; j < 3 * Nt; j++){
		current = Nodes[j];
		if (current != previous){
			counter++;
			previous = current;
			prAllNodes[counter] = previous;
		}
	}





	/*<<<<<<<<<<<<<<<<<<<<<<<< Assembly Phi matrix >>>>>>>>>>>>>>>>>>>>>>>>*/
	// Construct I, J, S for Phi matrix
	plhs[2] = mxCreateDoubleMatrix(2 * Phicount, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(2 * Phicount, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(2 * Phicount, 1, mxREAL);
	double* prI = mxGetPr(plhs[2]);
	double* prJ = mxGetPr(plhs[3]);
	double* prS = mxGetPr(plhs[4]);

	// Loop over all triangles and allocate barycentric coordinates for all atom sets
	counter = 0;
#pragma omp parallel
	{
		int _j, _m, _ind1, _ind2, _ind3;
		int* _prI = new int[2 * Phicount]; // private row indices
		int* _prJ = new int[2 * Phicount]; // private column indices
		double* _prS = new double[2 * Phicount]; // private data
		int _counter = 0; // private counter

#pragma omp for
		for (i = 0; i < Nt; i++){

			// Find Triangles[i].VertexAtoms[1:3] in prAllNodes
			_ind1 = -1;
			_ind2 = -1;
			_ind3 = -1;
			for (_m = 0; _m < Nrepatoms; _m++){
				if ((int)prAllNodes[_m] == (int)Triangles[i].VertexAtoms[1]){
					_ind1 = _m + 1; // MATLAB indexing
				}
				if ((int)prAllNodes[_m] == (int)Triangles[i].VertexAtoms[2]){
					_ind2 = _m + 1; // MATLAB indexing
				}
				if ((int)prAllNodes[_m] == (int)Triangles[i].VertexAtoms[3]){
					_ind3 = _m + 1; // MATLAB indexing
				}
				if ((_ind1 != -1) && (_ind2 != -1) && (_ind3 != -1)){
					break;
				}
			}
			if ((_ind1 == -1) || (_ind2 == -1) || (_ind3 == -1)){
				printf("repatoms array search failed.\n");
			}

			// Assign IntAtom's barycentric coordinates
			for (_j = 0; _j < Triangles[i].IntAtoms[0]; _j++){

				// Horizontal direction
				_prS[_counter + 0] = Triangles[i].Lint[3 * _j + 0];
				_prS[_counter + 1] = Triangles[i].Lint[3 * _j + 1];
				_prS[_counter + 2] = Triangles[i].Lint[3 * _j + 2];
				_prI[_counter + 0] = 2 * Triangles[i].IntAtoms[_j + 1] - 1;
				_prI[_counter + 1] = 2 * Triangles[i].IntAtoms[_j + 1] - 1;
				_prI[_counter + 2] = 2 * Triangles[i].IntAtoms[_j + 1] - 1;
				_prJ[_counter + 0] = 2 * _ind1 - 1;
				_prJ[_counter + 1] = 2 * _ind2 - 1;
				_prJ[_counter + 2] = 2 * _ind3 - 1;

				// Vertical direction
				_prS[_counter + 3] = Triangles[i].Lint[3 * _j + 0];
				_prS[_counter + 4] = Triangles[i].Lint[3 * _j + 1];
				_prS[_counter + 5] = Triangles[i].Lint[3 * _j + 2];
				_prI[_counter + 3] = 2 * Triangles[i].IntAtoms[_j + 1];
				_prI[_counter + 4] = 2 * Triangles[i].IntAtoms[_j + 1];
				_prI[_counter + 5] = 2 * Triangles[i].IntAtoms[_j + 1];
				_prJ[_counter + 3] = 2 * _ind1;
				_prJ[_counter + 4] = 2 * _ind2;
				_prJ[_counter + 5] = 2 * _ind3;

				_counter += 6;
			}

			// Assign EdgeAtom's barycentric coordinates
			for (_j = 0; _j < Triangles[i].EdgeAtoms[0]; _j++){

				// Horizontal direction
				_prS[_counter + 0] = Triangles[i].Ledge[3 * _j + 0];
				_prS[_counter + 1] = Triangles[i].Ledge[3 * _j + 1];
				_prS[_counter + 2] = Triangles[i].Ledge[3 * _j + 2];
				_prI[_counter + 0] = 2 * Triangles[i].EdgeAtoms[_j + 1] - 1;
				_prI[_counter + 1] = 2 * Triangles[i].EdgeAtoms[_j + 1] - 1;
				_prI[_counter + 2] = 2 * Triangles[i].EdgeAtoms[_j + 1] - 1;
				_prJ[_counter + 0] = 2 * _ind1 - 1;
				_prJ[_counter + 1] = 2 * _ind2 - 1;
				_prJ[_counter + 2] = 2 * _ind3 - 1;

				// Vertical direction
				_prS[_counter + 3] = Triangles[i].Ledge[3 * _j + 0];
				_prS[_counter + 4] = Triangles[i].Ledge[3 * _j + 1];
				_prS[_counter + 5] = Triangles[i].Ledge[3 * _j + 2];
				_prI[_counter + 3] = 2 * Triangles[i].EdgeAtoms[_j + 1];
				_prI[_counter + 4] = 2 * Triangles[i].EdgeAtoms[_j + 1];
				_prI[_counter + 5] = 2 * Triangles[i].EdgeAtoms[_j + 1];
				_prJ[_counter + 3] = 2 * _ind1;
				_prJ[_counter + 4] = 2 * _ind2;
				_prJ[_counter + 5] = 2 * _ind3;

				_counter += 6;
			}

			// Assign VertexAtom's barycentric coordinates
			if (3 == (int)Triangles[i].VertexAtoms[0]){

				// Horizontal direction
				_prS[_counter + 0] = 1;
				_prI[_counter + 0] = 2 * Triangles[i].VertexAtoms[1] - 1;
				_prJ[_counter + 0] = 2 * _ind1 - 1;
				_prS[_counter + 1] = 1;
				_prI[_counter + 1] = 2 * Triangles[i].VertexAtoms[2] - 1;
				_prJ[_counter + 1] = 2 * _ind2 - 1;
				_prS[_counter + 2] = 1;
				_prI[_counter + 2] = 2 * Triangles[i].VertexAtoms[3] - 1;
				_prJ[_counter + 2] = 2 * _ind3 - 1;

				// Vertical direction
				_prS[_counter + 3] = 1;
				_prI[_counter + 3] = 2 * Triangles[i].VertexAtoms[1];
				_prJ[_counter + 3] = 2 * _ind1;
				_prS[_counter + 4] = 1;
				_prI[_counter + 4] = 2 * Triangles[i].VertexAtoms[2];
				_prJ[_counter + 4] = 2 * _ind2;
				_prS[_counter + 5] = 1;
				_prI[_counter + 5] = 2 * Triangles[i].VertexAtoms[3];
				_prJ[_counter + 5] = 2 * _ind3;

				_counter += 6;
			}
		}

		// Collect _prI, _prJ, _prS from all threads
#pragma omp critical(ALLOCATE)
		{
			for (_j = 0; _j < _counter; _j++){
				prI[counter + _j] = _prI[_j];
				prJ[counter + _j] = _prJ[_j];
				prS[counter + _j] = _prS[_j];
			}
			counter += _counter;
			delete[] _prI;
			delete[] _prJ;
			delete[] _prS;
		}
	}

	// Resize outputs
	mxSetM(plhs[2], counter);
	mxSetM(plhs[3], counter);
	mxSetM(plhs[4], counter);

	// Delete Triangles database
	for (int j = 0; j < Nt; j++){
		delete[] Triangles[j].IntAtoms;
		delete[] Triangles[j].EdgeAtoms;
		delete[] Triangles[j].VertexAtoms;
		delete[] Triangles[j].Lint;
		delete[] Triangles[j].Ledge;
		delete[] Triangles[j].NeighTriangles;
	}
	delete[] Triangles;
	delete[] Nodes;
}

/******************** Function definitions ********************/
// Sort integers
int distInt(const void *a, const void *b){
	return(*(int*)a - *(int*)b);
}

// Compute centroid of a triangle according to flag
void TriagCentroid(double Cg[], double *P1, double *P2, double *P3, int flag){

	// Return centroid
	if (flag == 1){
		Cg[0] = (P1[0] + P2[0] + P3[0]) / 3;
		Cg[1] = (P1[1] + P2[1] + P3[1]) / 3;
	}

	// Return incenter
	else if (flag == 2){
		double a = sqrt(pow(P3[0] - P2[0], 2) + pow(P3[1] - P2[1], 2));
		double b = sqrt(pow(P3[0] - P1[0], 2) + pow(P3[1] - P1[1], 2));
		double c = sqrt(pow(P2[0] - P1[0], 2) + pow(P2[1] - P1[1], 2));
		Cg[0] = (a*P1[0] + b*P2[0] + c*P3[0]) / (a + b + c);
		Cg[1] = (a*P1[1] + b*P2[1] + c*P3[1]) / (a + b + c);
	}

	// Return circumcenter
	else if (flag == 3){
		double D = 2 * (P1[0] * (P2[1] - P3[1]) + P2[0] *
			(P3[1] - P1[1]) + P3[0] * (P1[1] - P2[1]));
		Cg[0] = ((pow(P1[0], 2) + pow(P1[1], 2))*(P2[1] - P3[1]) +
			(pow(P2[0], 2) + pow(P2[1], 2))*(P3[1] - P1[1]) +
			(pow(P3[0], 2) + pow(P3[1], 2))*(P1[1] - P2[1])) / D;
		Cg[1] = ((pow(P1[0], 2) + pow(P1[1], 2))*(P3[0] - P2[0]) +
			(pow(P2[0], 2) + pow(P2[1], 2))*(P1[0] - P3[0]) +
			(pow(P3[0], 2) + pow(P3[1], 2))*(P2[0] - P1[0])) / D;
	}
}
