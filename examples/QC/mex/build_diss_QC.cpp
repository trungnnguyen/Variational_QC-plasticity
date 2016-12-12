// build_diss_QC: computes the approximate dissipation part of the 
// incremental energy of the system

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

	// Test for five input arguments
	if (nrhs != 5){
		mexErrMsgTxt("Five input arguments required.");
	}

	// Get the data
	double* prz = mxGetPr(prhs[3]);
	double* przP = mxGetPr(prhs[4]);

	// Compute the dissipation distance
	nlhs = 1;
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* prDiss = mxGetPr(plhs[0]);
	int i;
	double Diss = 0;
	int NsAtoms = (int)mxGetNumberOfElements(prhs[1]);
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[2], "Potential");
	int name_ID = mxGetFieldNumber(prhs[1], "ID");
	int name_w = mxGetFieldNumber(prhs[1], "w");

#pragma omp parallel
	{
		int _j, _alpha, _bond;
		double *_BondList, *_prID, *_Potential, *_prw, _z_p, _Zp, _w;
		mxArray *_pBondList;

#pragma omp for reduction(+:Diss)
		for (i = 0; i < NsAtoms; i++){ // loop over all sampling atoms
			_prw = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_w));
			_prID = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_ID));
			_w = _prw[0];
			_alpha = (int)_prID[0]; // atom alpha ID, MATLAB indexing

			// Get data for an alpha atom
			_pBondList = mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_BondList);
			_BondList = mxGetPr(_pBondList);
			for (_j = 0; _j < mxGetN(_pBondList); _j++){ // loop over all the nearest neighbours of alpha sampling atom

				// Get data for the bond between atoms alpha and beta
				_bond = (int)_BondList[_j];
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[2],
					_bond - 1, name_Potential));
				_z_p = prz[_bond - 1]; // current plastic elongation
				_Zp = przP[_bond - 1]; // previous time-step plastic elongation

				// Compute D part of the energy
				Diss += _w*(0.5 * _Potential[2] * fabs(_z_p - _Zp));
			}
		}
	}

	prDiss[0] = Diss;
}
