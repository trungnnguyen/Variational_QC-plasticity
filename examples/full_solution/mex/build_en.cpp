// build_en: computes the elastic part of the incremental energy of the
// system

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

	// Test for five input argument
	if (nrhs != 5){
		mexErrMsgTxt("Five input arguments required.");
	}

	// Get the input data
	double* prz = mxGetPr(prhs[2]);
	double* prr = mxGetPr(prhs[3]);
	double* prk = mxGetPr(prhs[4]);

	// Compute energy
	nlhs = 1;
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* prEn = mxGetPr(plhs[0]);
	int i;
	double En = 0.0;
	int NAtoms = (int)mxGetNumberOfElements(prhs[0]);
	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[1], "Potential");

#pragma omp parallel
	{
		int _j, _alpha, _beta, _bond;
		double _R, _r, _EnInt, _EnHard, _z_p, _z_c;
		double *_Ra, *_Rb, *_NeighbourList, *_BondList, *_Potential;
		double _ra[2], _rb[2];
		mxArray *_pNeighbourList; // pointer to NeighbourList array

#pragma omp for reduction(+:En)
		for (i = 0; i < NAtoms; i++){ // loop over all atoms
			_alpha = i + 1; // MATLAB indexing

			// Get data for atom alpha
			_Ra = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1, name_R));
			_pNeighbourList = mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_NeighbourList);
			_NeighbourList = mxGetPr(_pNeighbourList);
			_BondList = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_BondList));
			_ra[0] = prr[2 * (_alpha - 1)];
			_ra[1] = prr[2 * (_alpha - 1) + 1];
			for (_j = 0; _j < mxGetN(_pNeighbourList); _j++){ // loop over all nearest neighbours to atom alpha
				_beta = (int)_NeighbourList[_j]; // atom beta ID, MATLAB indexing

				// Get data for beta atom
				_Rb = mxGetPr(mxGetFieldByNumber(prhs[0], _beta - 1,
					name_R));
				_rb[0] = prr[2 * (_beta - 1)];
				_rb[1] = prr[2 * (_beta - 1) + 1];
				_R = sqrt(pow(_Rb[0] - _Ra[0], 2) +
					pow(_Rb[1] - _Ra[1], 2)); // the initial bond lengtht
				_r = sqrt(pow(_rb[0] - _ra[0], 2) +
					pow(_rb[1] - _ra[1], 2)); // bond length

				// Get data for the bond between atoms alpha and beta
				_bond = (int)_BondList[_j];
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[1], _bond - 1,
					name_Potential));
				_z_p = prz[_bond - 1]; // current plastic elongation
				_z_c = prk[_bond - 1]; // current cumulative plastic elongation

				// Compute internal part of the energy
				_EnInt = 0.25 * (_Potential[0] / _R)*pow(_r - _R - _z_p, 2);

				// Compute hardening part of the energy
				_EnHard = 0.5 * (1 / (_Potential[3] + 1)) * _Potential[1] *
					_Potential[2] * _R * pow(_z_c / _R, _Potential[3] + 1);

				// Add to the overall energy
				En += (_EnInt + _EnHard);
			}
		}
	}

	// Send the data back to MATLAB
	prEn[0] = En;
}
