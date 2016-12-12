// sort_sampling_atoms_QC: determine sampling atoms and their weights.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "mex.h"
#include "matrix.h"

/******************** Structure declarations ******************/
struct SATOM{
	int ID; // ID of sampling atom in atoms database
	int Nneigh; // number of neighbours within the same triangle
	double w; // sampling atom weight, possibly non-integer
	int* List; // list of represented atoms
	double r; // distance from triangle's centre - for selecting an appropriate satom
};

/******************** Function declarations *******************/
int distAtoms(const void *AtomA, const void *AtomB); // comparison function for distance of atoms
int distInt(const void *a, const void *b); // comparison function for integers
void TriagCentroid(double Cg[], double *P1, double *P2, double *P3, int flag); // according to flag compute internal point of a triangle
int distAtomsNneigh(const void *a, const void *b); // comparison function for distance of atoms according to Nneigh
int distAtomsR(const void *a, const void *b); // comparison function for distance of atoms according to r
int LengthIntersect(int Na, int *A, int Nb, int *B); // for sorted arrays A, B returns length(intersect(A,B))

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for six input arguments
	if (nrhs != 6){
		mexErrMsgTxt("Six input arguments required.");
	}

	// Copy input variables
	int NAtoms = (int)mxGetNumberOfElements(prhs[0]);
	int NTriangles = (int)mxGetNumberOfElements(prhs[1]);
	double* prSumRule = mxGetPr(prhs[2]);
	int SumRule = (int)prSumRule[0];
	double* prrhs = mxGetPr(prhs[3]);
	double SizeX = prrhs[0];
	double SizeY = prrhs[1];
	double *prSW = mxGetPr(prhs[4]);
	int SW = (int)prSW[0];
	double* prTOL = mxGetPr(prhs[5]);
	double TOL = prTOL[0]; // a positive number < TOL is treated as zero         

	// Allocate memory
	int i, j;
	int maxS = NAtoms;
	SATOM* sampatoms = new SATOM[maxS];
	if (sampatoms == NULL){
		mexErrMsgTxt("Not enough memory for sampling atoms.");
	}

	/*<<<<<<<<<<<<<<<<<<<<<<<< Full summation rule >>>>>>>>>>>>>>>>>>>>>>>>*/
	// Add all atoms
	int countersampatoms;
	if (SumRule == 2){
		for (i = 0; i < NAtoms; i++){
			sampatoms[i].ID = i + 1; // MATLAB indexing
			sampatoms[i].w = 1;
			sampatoms[i].r = 0;
			sampatoms[i].List = new int[2];
			sampatoms[i].List[0] = 1;
			sampatoms[i].List[1] = i + 1;
		}
		countersampatoms = NAtoms;
	}




	/*<<<<<<<<<<<<<<<<<<<<<<<< Exact or central summation rule >>>>>>>>>>>>>>>>>>>>>>>>*/
	else if ((SumRule == 0) || (SumRule == 1)){

		int maxIntAtoms = 0, maxEdgeAtoms = 0;
		double *IntAtoms, *EdgeAtoms;
		mxArray *pIntAtoms, *pEdgeAtoms;
		int name_IntAtoms = mxGetFieldNumber(prhs[1], "IntAtoms");
		int name_EdgeAtoms = mxGetFieldNumber(prhs[1], "EdgeAtoms");

		// Find maximum number of IntAtoms within a triangle
		for (i = 0; i < NTriangles; i++){
			pIntAtoms = mxGetFieldByNumber(prhs[1], i, name_IntAtoms);
			IntAtoms = mxGetPr(pIntAtoms);
			if (maxIntAtoms < mxGetN(pIntAtoms)){
				maxIntAtoms = (int)mxGetN(pIntAtoms);
			}
		}

		// Find maximum number of EdgeAtoms within a triangle
		for (i = 0; i < NTriangles; i++){
			pEdgeAtoms = mxGetFieldByNumber(prhs[1], i, name_EdgeAtoms);
			EdgeAtoms = mxGetPr(pEdgeAtoms);
			if (maxEdgeAtoms < mxGetN(pEdgeAtoms)){
				maxEdgeAtoms = (int)mxGetN(pEdgeAtoms);
			}
		}

		// Allocate variables
		int m, n, k, counterList, alpha, tNeighbourList[8], tVertexAtoms[3];
		int nNeighbourList, nVertexAtoms, nNeighTriangles, ntempIntAtoms,
			ntempEdgeAtoms, nIntAtoms, nEdgeAtoms;
		int* List = new int[maxIntAtoms + maxEdgeAtoms];
		double *Ra, *T, we;
		double *NeighbourList, *VertexAtoms, *NeighTriangles, *tempIntAtoms,
			*tempEdgeAtoms;
		SATOM* IntData = new SATOM[maxIntAtoms]; // working database of internal atoms
		mxArray *pNeighbourList, *pVertexAtoms, *pNeighTriangles, *ptempIntAtoms,
			*ptempEdgeAtoms;
		int *tIntAtoms = new int[maxIntAtoms];
		int *tEdgeAtoms = new int[maxEdgeAtoms];
		int *ttempEdgeAtoms = new int[maxEdgeAtoms];
		int name_T = mxGetFieldNumber(prhs[1], "T");
		int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
		int name_VertexAtoms = mxGetFieldNumber(prhs[1], "VertexAtoms");
		int name_NeighTriangles = mxGetFieldNumber(prhs[1], "NeighTriangles");
		int name_R = mxGetFieldNumber(prhs[0], "R");
		countersampatoms = 0;

		// For every triangle find its sampling atoms
		for (i = 0; i < NTriangles; i++){
			counterList = 0;
			pIntAtoms = mxGetFieldByNumber(prhs[1], i, name_IntAtoms);
			IntAtoms = mxGetPr(pIntAtoms);
			nIntAtoms = (int)mxGetN(pIntAtoms);
			pEdgeAtoms = mxGetFieldByNumber(prhs[1], i, name_EdgeAtoms);
			EdgeAtoms = mxGetPr(pEdgeAtoms);
			nEdgeAtoms = (int)mxGetN(pEdgeAtoms);
			pVertexAtoms = mxGetFieldByNumber(prhs[1], i, name_VertexAtoms);
			VertexAtoms = mxGetPr(pVertexAtoms);
			nVertexAtoms = (int)mxGetN(pVertexAtoms);
			pNeighTriangles = mxGetFieldByNumber(prhs[1], i,
				name_NeighTriangles);
			NeighTriangles = mxGetPr(pNeighTriangles);
			nNeighTriangles = (int)mxGetN(pNeighTriangles);

			// Add all three VertexAtoms to sampatoms
			for (j = 0; j < nVertexAtoms; j++){
				sampatoms[countersampatoms].ID = (int)VertexAtoms[j]; // MATLAB indexing
				sampatoms[countersampatoms].w = 1;
				sampatoms[countersampatoms].r = 0;
				sampatoms[countersampatoms].List = new int[2];
				sampatoms[countersampatoms].List[0] = 1;
				sampatoms[countersampatoms].List[1] = (int)VertexAtoms[j];
				countersampatoms++;
			}

			// Sort nIntAtoms according to the distance from T
			if (nIntAtoms > 0){
				T = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_T));
				for (k = 0; k < nIntAtoms; k++){
					alpha = (int)IntAtoms[k];
					Ra = mxGetPr(mxGetFieldByNumber(prhs[0],
						alpha - 1, name_R));
					IntData[k].ID = alpha; // MATLAB indexing
					IntData[k].r = sqrt(pow(Ra[0] - T[0], 2) +
						pow(Ra[1] - T[1], 2)); // distance of the satom from the triangle incenter

				}
				qsort(IntData, k, sizeof(SATOM), distAtomsR); // sort by distance
			}





			/*<<<<<<<<<<<<<<<<<<<<<<<< Exact summation rule >>>>>>>>>>>>>>>>>>>>>>>>*/
			if (SumRule == 0){

				// Add all EdgeAtoms to sampatoms
				for (j = 0; j < nEdgeAtoms; j++){
					sampatoms[countersampatoms].ID = (int)EdgeAtoms[j]; // MATLAB indexing
					sampatoms[countersampatoms].w = 1;
					sampatoms[countersampatoms].r = 0;
					sampatoms[countersampatoms].List = new int[2];
					sampatoms[countersampatoms].List[0] = 1;
					sampatoms[countersampatoms].List[1] = (int)EdgeAtoms[j];
					countersampatoms++;
				}

				if (nIntAtoms > 0){

					// Assign values to temporary vectors
					for (j = 0; j < nIntAtoms; j++){
						tIntAtoms[j] = (int)IntAtoms[j];
					}
					for (j = 0; j < nVertexAtoms; j++){
						tVertexAtoms[j] = (int)VertexAtoms[j];
					}
					for (j = 0; j < nEdgeAtoms; j++){
						tEdgeAtoms[j] = (int)EdgeAtoms[j];
					}

					// Sort temporary vectors
					qsort(tIntAtoms, nIntAtoms, sizeof(int), distInt);
					qsort(tVertexAtoms, nVertexAtoms, sizeof(int), distInt);
					qsort(tEdgeAtoms, nEdgeAtoms, sizeof(int), distInt);

					// Find those atoms from IntAtoms that do not have all the nearest neighbours within the current triangle
					for (j = 0; j < nIntAtoms; j++){
						alpha = IntData[j].ID; // MATLAB indexing
						pNeighbourList = mxGetFieldByNumber(prhs[0],
							alpha - 1, name_NeighbourList);
						NeighbourList = mxGetPr(pNeighbourList); // get neighbourlist of alpha
						nNeighbourList = (int)mxGetN(pNeighbourList);

						// Temporary NeighbourList
						for (m = 0; m < nNeighbourList; m++){
							tNeighbourList[m] = (int)NeighbourList[m];
						}
						qsort(tNeighbourList, nNeighbourList, sizeof(int),
							distInt);
						IntData[j].Nneigh = LengthIntersect(nIntAtoms,
							tIntAtoms, nNeighbourList, tNeighbourList); // length(intersect(tIntAtoms,tNeighbourList))
						IntData[j].Nneigh += LengthIntersect(nEdgeAtoms,
							tEdgeAtoms, nNeighbourList, tNeighbourList); // length(intersect(tEdgeAtoms,tNeighbourList))
						IntData[j].Nneigh += LengthIntersect(nVertexAtoms,
							tVertexAtoms, nNeighbourList, tNeighbourList); // length(intersect(tVertexAtoms,tNeighbourList))

						// If not all neighbours are inside the triangle, add this atom into sampatoms list as a discrete sampling atom
						if (IntData[j].Nneigh < nNeighbourList){
							sampatoms[countersampatoms].ID = alpha; // MATLAB indexing
							sampatoms[countersampatoms].w = 1;
							sampatoms[countersampatoms].r = 0;
							sampatoms[countersampatoms].List = new int[2];
							sampatoms[countersampatoms].List[0] = 1;
							sampatoms[countersampatoms].List[1] = alpha;
							countersampatoms++;
						}

						// If all the neighbours are inside the triangle, add this atom to List of the central sampling atom
						else{
							List[counterList] = alpha; // List of atoms with all the neighbours within a triangle
							counterList++;
						}
					}

					// Update central sampling atom database to include the central sampling atom
					if (counterList != 0){
						sampatoms[countersampatoms].ID = List[0]; // MATLAB indexing
						sampatoms[countersampatoms].w = counterList;
						sampatoms[countersampatoms].r = 0;
						sampatoms[countersampatoms].List =
							new int[counterList + 1];
						sampatoms[countersampatoms].List[0] = counterList;

						// Copy List into the corresponding sampatoms List
						for (n = 0; n < counterList; n++){
							sampatoms[countersampatoms].List[n + 1] = List[n]; // List of atoms with all the neighbours within a triangle
						}
						countersampatoms++;
					}
				}
			}





			/*<<<<<<<<<<<<<<<<<<<<<<<< Central summation rule >>>>>>>>>>>>>>>>>>>>>>>>*/
			else if (SumRule == 1){

				// If there is no internal atom inside the triangle, then add all edge atoms as discrete sampling atoms
				if (nIntAtoms == 0){

					// Add all EdgeAtoms to sampatoms
					for (j = 0; j < nEdgeAtoms; j++){
						alpha = (int)EdgeAtoms[j]; // MATLAB indexing
						sampatoms[countersampatoms].ID = alpha;
						sampatoms[countersampatoms].w = 1;
						sampatoms[countersampatoms].r = 0;
						sampatoms[countersampatoms].List = new int[2];
						sampatoms[countersampatoms].List[0] = 1;
						sampatoms[countersampatoms].List[1] = alpha;
						countersampatoms++;
					}
				}

				// If there is at least one internal atom, choose the central sampling atom and compute its weight
				else{

					// Check the number of the nearest neighbours within a triangle for the atom closest to T
					alpha = IntData[0].ID; // MATLAB indexing
					pNeighbourList = mxGetFieldByNumber(prhs[0],
						alpha - 1, name_NeighbourList);
					NeighbourList = mxGetPr(pNeighbourList); // get the neighbourlist of alpha
					nNeighbourList = (int)mxGetN(pNeighbourList);

					// Assign values to temporary vectors
					for (j = 0; j < nIntAtoms; j++){
						tIntAtoms[j] = (int)IntAtoms[j];
					}
					for (j = 0; j < nVertexAtoms; j++){
						tVertexAtoms[j] = (int)VertexAtoms[j];
					}
					for (j = 0; j < nEdgeAtoms; j++){
						tEdgeAtoms[j] = (int)EdgeAtoms[j];
					}
					for (j = 0; j < nNeighbourList; j++){
						tNeighbourList[j] = (int)NeighbourList[j];
					}

					// Sort the temporary vectors
					qsort(tIntAtoms, nIntAtoms, sizeof(int), distInt);
					qsort(tVertexAtoms, nVertexAtoms, sizeof(int), distInt);
					qsort(tEdgeAtoms, nEdgeAtoms, sizeof(int), distInt);
					qsort(tNeighbourList, nNeighbourList, sizeof(int),
						distInt);

					// Number of the nearest neighbours inside triangle
					IntData[0].Nneigh = LengthIntersect(nIntAtoms,
						tIntAtoms, nNeighbourList, tNeighbourList); // length(intersect(tIntAtoms,tNeighbourList))
					IntData[0].Nneigh += LengthIntersect(nEdgeAtoms,
						tEdgeAtoms, nNeighbourList, tNeighbourList); // length(intersect(tEdgeAtoms,tNeighbourList))
					IntData[0].Nneigh += LengthIntersect(nVertexAtoms,
						tVertexAtoms, nNeighbourList, tNeighbourList); // length(intersect(tVertexAtoms,tNeighbourList))

					// If atom alpha (closest to T) has less neighbours within the triangle than 8, or if it is not the only internal atom of the triangle, then verify all the remaining IntAtoms
					if ((IntData[0].Nneigh != 8) || (nIntAtoms != 1)){

						// For all IntAtoms compute the number of neighbours within the triangle, Nneigh
						for (k = 1; k < nIntAtoms; k++){
							alpha = IntData[k].ID;
							pNeighbourList = mxGetFieldByNumber(prhs[0],
								alpha - 1, name_NeighbourList);
							NeighbourList = mxGetPr(pNeighbourList); // get the neighbourlist of alpha
							nNeighbourList = (int)mxGetN(pNeighbourList);

							// Compute the number of neighbours inside triangle, Nneigh
							for (j = 0; j < nNeighbourList; j++){
								tNeighbourList[j] = (int)NeighbourList[j];
							}
							qsort(tNeighbourList, nNeighbourList,
								sizeof(int), distInt);
							IntData[k].Nneigh = LengthIntersect(nIntAtoms,
								tIntAtoms, nNeighbourList, tNeighbourList); // length(intersect(tIntAtoms,tNeighbourList))
							IntData[k].Nneigh += LengthIntersect(nEdgeAtoms,
								tEdgeAtoms, nNeighbourList, tNeighbourList); // length(intersect(tEdgeAtoms,tNeighbourList))
							IntData[k].Nneigh += LengthIntersect(nVertexAtoms,
								tVertexAtoms, nNeighbourList, tNeighbourList); // length(intersect(tVertexAtoms,tNeighbourList))
						}

						// Sort IntData according to number of neighbours inside triangle, Nneigh
						qsort(IntData, nIntAtoms, sizeof(SATOM),
							distAtomsNneigh);

						// Count the number of competitors having the same number of neighbours inside the triangle
						k = 1;
						while (IntData[0].Nneigh == IntData[k].Nneigh){
							k++;
						}

						// Sort the first k competitors according to r
						qsort(IntData, k, sizeof(SATOM), distAtomsR);
					}

					// Assembly corresponding weight and List
					counterList = 0;
					we = 0;

					// For each neighbouring triangle with shared edge
					for (j = 0; j < nNeighTriangles; j++){
						ptempIntAtoms = mxGetFieldByNumber(prhs[1],
							(int)NeighTriangles[j] - 1, name_IntAtoms);
						tempIntAtoms = mxGetPr(ptempIntAtoms);
						ntempIntAtoms = (int)mxGetN(ptempIntAtoms);

						// If neighbouring triangle has an internal atom, count the number of atoms lying on a common edge
						if (ntempIntAtoms > 0){

							// Find intersect(EdgeAtoms,ttempEdgeAtoms) and assign it to List
							ptempEdgeAtoms = mxGetFieldByNumber(prhs[1],
								(int)NeighTriangles[j] - 1, name_EdgeAtoms);
							tempEdgeAtoms = mxGetPr(ptempEdgeAtoms);
							ntempEdgeAtoms = (int)mxGetN(ptempEdgeAtoms);
							for (m = 0; m < ntempEdgeAtoms; m++){
								ttempEdgeAtoms[m] = (int)tempEdgeAtoms[m];
							}
							qsort(ttempEdgeAtoms, ntempEdgeAtoms,
								sizeof(int), distInt);

							// Add intersection(ttempEdgeAtoms,tEdgeAtoms) to the List
							m = 0;
							n = 0;
							while (m != ntempEdgeAtoms && n != nEdgeAtoms)
							{
								if (ttempEdgeAtoms[m] < tEdgeAtoms[n]) m++;
								else if (tEdgeAtoms[n] < ttempEdgeAtoms[m]) n++;
								else {
									we += 0.5;
									List[counterList] = ttempEdgeAtoms[m];
									counterList++;
									m++;
									n++;
								}
							}
						}

						// Else do nothing, since atoms lying on the common
						// edge with a neighbouring triangle that has no 
						// internal atom are discrete sampling atoms
					}

					// Add atoms lying on the physical boundary
					for (m = 0; m < nEdgeAtoms; m++){
						alpha = (int)EdgeAtoms[m]; // MATLAB indexing
						Ra = mxGetPr(mxGetFieldByNumber(prhs[0],
							alpha - 1, name_R));

						// If alpha lies on the physical boundary, add it to the List -- possibly choose another value than 0.5 for its weight
						if (SW == 0){ // uniform loading test
							if ((fabs(fabs(Ra[0]) - SizeX) < TOL) ||
								(fabs(fabs(Ra[1]) - SizeY) < TOL)){ // holds for the entire boundary
								we += 0.5;
								List[counterList] = alpha;
								counterList++;
							}
						}
						else if (SW == 1){ // pure bending test
							if (fabs(fabs(Ra[1]) - SizeY) < TOL){ // the case of the Gamma_1 and Gamma_3 parts of the boundary
								we += 0.5;
								List[counterList] = alpha;
								counterList++;
							}
							if (fabs(fabs(Ra[0]) - SizeX) < TOL){ // the case of Gamma_2 and Gamma_4 parts of the boundary
								we += 0.5;
								List[counterList] = alpha;
								counterList++;
							}
						}
					}

					// Copy extracted data to sampatoms database
					// From k-competitors, add to sampatoms the closest one to the triangle's internal point T defined in sort_atoms_QC
					sampatoms[countersampatoms].ID = IntData[0].ID; // MATLAB indexing
					sampatoms[countersampatoms].w = nIntAtoms + we;
					sampatoms[countersampatoms].r = 0;
					sampatoms[countersampatoms].List = new int[nIntAtoms +
						counterList + 1];
					sampatoms[countersampatoms].List[0] = nIntAtoms +
						counterList;

					// Copy all the IntAtoms into the corresponding sampling atom
					for (n = 0; n < nIntAtoms; n++){
						sampatoms[countersampatoms].List[n + 1] =
							(int)IntAtoms[n];
					}

					// Copy all the appropriate EdgeAtoms from List into the corresponding sampling atom databases
					for (n = 0; n < counterList; n++){
						sampatoms[countersampatoms].List[nIntAtoms + n + 1] =
							List[n];
					}
					countersampatoms++;
				}
			}





			/*<<<<<<<<<<<<<<<<<<<<<<<< Reallocate sampatoms database, if necessary >>>>>>>>>>>>>>>>>>>>>>>>*/
			if (countersampatoms > (int)(0.8*maxS)){
				maxS += NAtoms;
				SATOM* temp1 = new SATOM[maxS];
				for (n = 0; n < countersampatoms; n++){
					temp1[n].ID = sampatoms[n].ID;
					temp1[n].w = sampatoms[n].w;
					temp1[n].r = sampatoms[n].r;
					temp1[n].List = new int[sampatoms[n].List[0] + 1];
					for (k = 0; k < sampatoms[n].List[0] + 1; k++){
						temp1[n].List[k] = sampatoms[n].List[k];
					}
				}

				// Delete sampatoms
				for (n = 0; n < countersampatoms; n++){
					delete[] sampatoms[n].List;
				}
				delete[] sampatoms;
				sampatoms = temp1;
			}
		}
		delete[] List;
		delete[] IntData;
		delete[] tIntAtoms;
		delete[] tEdgeAtoms;
		delete[] ttempEdgeAtoms;
	}





	/*<<<<<<<<<<<<<<<<<<<<<<<< Send out the data back to MATLAB >>>>>>>>>>>>>>>>>>>>>>>>*/
	// Sort sampatoms according to their IDs, remove zeros
	qsort(sampatoms, countersampatoms, sizeof(SATOM), distAtoms);

	// Count size of the output
	int previous = sampatoms[0].ID;
	int MOut = 0, current;
	for (i = 0; i < countersampatoms; i++){ // test for multiple ID occurences in sampatoms database
		current = sampatoms[i].ID;
		if (current != previous){
			previous = current;
			MOut++;
		}
	}
	nlhs = 1;
	mwSize dims[] = { 1, MOut + 1 };
	const char *field_names[] = { "ID", "w", "List" };
	plhs[0] = mxCreateStructArray(2, dims, 3, field_names);

	// Populate output with computed data
	int name_ID = mxGetFieldNumber(plhs[0], "ID");
	int name_w = mxGetFieldNumber(plhs[0], "w");
	int name_List = mxGetFieldNumber(plhs[0], "List");
	previous = sampatoms[0].ID;
	int counter = 0;
	for (i = 0; i < countersampatoms; i++){
		current = sampatoms[i].ID;

		// Test for multiple ID occurences in sampatoms database
		if (current != previous || i == 0){
			mxArray *array_ID, *array_w, *array_List;
			double *prID, *prw, *prList;

			// Create arrays
			array_ID = mxCreateDoubleMatrix(1, 1, mxREAL);
			array_w = mxCreateDoubleMatrix(1, 1, mxREAL);
			array_List = mxCreateDoubleMatrix(1, sampatoms[i].List[0],
				mxREAL);

			// Get data arrays
			prID = mxGetPr(array_ID);
			prw = mxGetPr(array_w);
			prList = mxGetPr(array_List);

			// Populate data arrays
			prID[0] = sampatoms[i].ID;
			prw[0] = sampatoms[i].w;
			for (j = 0; j < sampatoms[i].List[0]; j++){
				prList[j] = (double)sampatoms[i].List[j + 1];
			}

			// Assign arrays
			mxSetFieldByNumber(plhs[0], counter, name_ID, array_ID);
			mxSetFieldByNumber(plhs[0], counter, name_w, array_w);
			mxSetFieldByNumber(plhs[0], counter, name_List, array_List);
			previous = current;
			counter++;
		}
	}

	// Delete sampatoms
	for (i = 0; i < countersampatoms; i++){
		delete[] sampatoms[i].List;
	}
	delete[] sampatoms;
}

/******************** Function definitions ********************/
// Sort integers
int distInt(const void *a, const void *b){
	return(*(int*)a - *(int*)b);
}

// Sort atoms with respect to ID
int distAtoms(const void *a, const void *b){
	SATOM *A = (SATOM*)a,
		*B = (SATOM*)b;
	return(A->ID - B->ID);
}

// Sort atoms with respect to weights
int distAtomsNneigh(const void *a, const void *b){
	SATOM *A = (SATOM*)a,
		*B = (SATOM*)b;
	if (A->Nneigh > B->Nneigh) return -1;
	else if (A->Nneigh < B->Nneigh) return 1;
	else return 0;
}

// Sort atoms with respect to distance
int distAtomsR(const void *a, const void *b){
	SATOM *A = (SATOM*)a,
		*B = (SATOM*)b;
	if (A->r < B->r) return -1;
	else if (A->r > B->r) return 1;
	else return 0;
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
		double D = 2 * (P1[0] * (P2[1] - P3[1]) + P2[0] * (P3[1] - P1[1])
			+ P3[0] * (P1[1] - P2[1]));
		Cg[0] = ((pow(P1[0], 2) + pow(P1[1], 2))*(P2[1] - P3[1]) +
			(pow(P2[0], 2) + pow(P2[1], 2))*(P3[1] - P1[1]) +
			(pow(P3[0], 2) + pow(P3[1], 2))*(P1[1] - P2[1])) / D;
		Cg[1] = ((pow(P1[0], 2) + pow(P1[1], 2))*(P3[0] - P2[0]) +
			(pow(P2[0], 2) + pow(P2[1], 2))*(P1[0] - P3[0]) +
			(pow(P3[0], 2) + pow(P3[1], 2))*(P2[0] - P1[0])) / D;
	}
}

// Compute length(intersect(A,B))
int LengthIntersect(int Na, int *A, int Nb, int *B){
	int i = 0;
	int j = 0;
	int count = 0;
	while (i != Na && j != Nb)
	{
		if (A[i] < B[j]) i++;
		else if (B[j] < A[i]) j++;
		else {
			count++;
			i++;
			j++;
		}
	}
	return count;
}
