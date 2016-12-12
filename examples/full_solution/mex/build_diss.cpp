// build_diss: computes the dissipation part of the incremental energy of
// the system

#include <stdio.h>
#include <math.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mex.h"
#include "matrix.h"

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for four input argument
	if (nrhs != 4){
		mexErrMsgTxt("Four input arguments required.");
	}

	// Get the input data
	double* prz = mxGetPr(prhs[2]);
	double* prk = mxGetPr(prhs[3]);

	// Compute the dissipation distance
	nlhs = 1;
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* prDiss = mxGetPr(plhs[0]);
	int i;
	double Diss = 0;
	int NAtoms = (int)mxGetNumberOfElements(prhs[0]);
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[1], "Potential");

#pragma omp parallel
	{
		int _j, _alpha, _bond;
		double *_BondList, *_Potential, _z_p, _Zp;
		mxArray *_pBondList; // pointer to BondList array

#pragma omp for reduction(+:Diss)
		for (i = 0; i < NAtoms; i++){ // loop over all atoms
			_alpha = i + 1; // MATLAB indexing

			// Get data for atom alpha
			_pBondList = mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_BondList);
			_BondList = mxGetPr(_pBondList);
			for (_j = 0; _j < mxGetN(_pBondList); _j++){ // loop over all nearest neighbours to atom alpha

				// Get data for the bond between atoms alpha and beta
				_bond = (int)_BondList[_j];
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[1], _bond - 1,
					name_Potential));
				_z_p = prz[_bond - 1]; // current plastic elongation
				_Zp = prk[_bond - 1]; // previous-time-step plastic elongation

				// Compute dissipation distance D for a single bond
				Diss += (0.5 * _Potential[2] * fabs(_z_p - _Zp));
			}
		}
	}

	prDiss[0] = Diss;
}
