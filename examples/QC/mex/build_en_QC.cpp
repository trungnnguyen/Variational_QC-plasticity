// build_en_QC: computes the approximate elastic part of the incremental
// energy of the system

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

	// Test for six input arguments
	if (nrhs != 6){
		mexErrMsgTxt("Six input arguments required.");
	}

	// Get the data
	double* prz = mxGetPr(prhs[3]);
	double* prr = mxGetPr(prhs[4]);
	double* prk = mxGetPr(prhs[5]);

	// Compute energy
	nlhs = 1;
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* prEn = mxGetPr(plhs[0]);
	int i;
	double En = 0.0;
	int NsAtoms = (int)mxGetNumberOfElements(prhs[1]);
	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[2], "Potential");
	int name_ID = mxGetFieldNumber(prhs[1], "ID");
	int name_w = mxGetFieldNumber(prhs[1], "w");

#pragma omp parallel
	{
		int _j, _alpha, _beta, _bond;
		double *_Ra, *_Rb, *_NeighbourList, *_BondList, *_prID,
			*_Potential, *_prw;
		double _ra[2], _rb[2], _R, _r, _z_p, _z_c, _w, _EnInt, _EnHard;
		mxArray *_pNeighbourList;

#pragma omp for reduction(+:En)
		for (i = 0; i < NsAtoms; i++){ // loop over all sampling atoms
			_prw = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_w));
			_prID = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_ID));
			_w = _prw[0];
			_alpha = (int)_prID[0]; // atom alpha ID, MATLAB indexing

			// Get data for an alpha sampling atom
			_Ra = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1, name_R));
			_pNeighbourList = mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_NeighbourList);
			_NeighbourList = mxGetPr(_pNeighbourList);
			_BondList = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_BondList));
			_ra[0] = prr[2 * (_alpha - 1)];
			_ra[1] = prr[2 * (_alpha - 1) + 1];
			for (_j = 0; _j < mxGetN(_pNeighbourList); _j++){ // loop over all nearest neighbours of sampling atom alpha
				_beta = (int)_NeighbourList[_j]; // atom beta ID, MATLAB indexing

				// Get data for beta atom
				_Rb = mxGetPr(mxGetFieldByNumber(prhs[0], _beta - 1,
					name_R));
				_rb[0] = prr[2 * (_beta - 1)];
				_rb[1] = prr[2 * (_beta - 1) + 1];
				_R = sqrt(pow(_Rb[0] - _Ra[0], 2) +
					pow(_Rb[1] - _Ra[1], 2)); // the original bond length
				_r = sqrt(pow(_rb[0] - _ra[0], 2) +
					pow(_rb[1] - _ra[1], 2)); // bond length

				// Get data for the bond between atoms alpha and beta
				_bond = (int)_BondList[_j];
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[2], _bond - 1,
					name_Potential));
				_z_p = prz[_bond - 1]; // current plastic elongation
				_z_c = prk[_bond - 1]; // current cumulative plastic elongation

				// Compute internal part of the energy
				_EnInt = 0.25 * (_Potential[0] / _R)*pow(_r - _R - _z_p, 2);

				// Compute hardening part of the energy
				_EnHard = 0.5 * (1 / (_Potential[3] + 1)) * _Potential[1] *
					_Potential[2] * _R * pow(_z_c / _R, _Potential[3] + 1);

				// Add to overall energy
				En += _w*(_EnInt + _EnHard);
			}
		}
	}

	// Send the data back to MATLAB
	prEn[0] = En;
}
