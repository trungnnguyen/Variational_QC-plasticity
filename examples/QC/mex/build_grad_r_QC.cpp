// build_grad_r_QC: assembly the gradient of approximate internal
// energy w.r.t. r using summation rule

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

	// Get the input data
	double* prr = mxGetPr(prhs[3]);
	double* prz = mxGetPr(prhs[4]);
	int NAtoms = (int)mxGetNumberOfElements(prhs[0]);
	int NsAtoms = (int)mxGetNumberOfElements(prhs[1]);

	// Compute the gradient and send out the data back to MATLAB
	nlhs = 1;
	mxArray* Grad = mxCreateDoubleMatrix(2 * NAtoms, 1, mxREAL);
	double* prGrad = mxGetPr(Grad);
	plhs[0] = Grad;
	int i;
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
		double _ra[2], _rb[2], _rab[2], _fa[2], _R, _r, _dphi, _z_p, _w;
		mxArray *_pNeighbourList;

#pragma omp for
		for (i = 0; i < NsAtoms; i++){ // loop over all sampling atoms
			_prw = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_w));
			_prID = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_ID));
			_w = _prw[0];
			_alpha = (int)_prID[0]; // atom alpha ID, MATLAB indexing

			// Get data for an i-th sampling atom, atom alpha
			_Ra = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1, name_R));
			_pNeighbourList = mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_NeighbourList);
			_NeighbourList = mxGetPr(_pNeighbourList);
			_BondList = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_BondList));
			_ra[0] = prr[2 * (_alpha - 1)];
			_ra[1] = prr[2 * (_alpha - 1) + 1];
			for (_j = 0; _j < mxGetN(_pNeighbourList); _j++){ // loop over all nearest neighbours
				_beta = (int)_NeighbourList[_j]; // atom beta ID, MATLAB indexing

				// Get data for j-th neighboring atom, atom beta
				_Rb = mxGetPr(mxGetFieldByNumber(prhs[0], _beta - 1,
					name_R));
				_rb[0] = prr[2 * (_beta - 1)];
				_rb[1] = prr[2 * (_beta - 1) + 1];
				_R = sqrt(pow(_Rb[0] - _Ra[0], 2) +
					pow(_Rb[1] - _Ra[1], 2)); // the initial bond length
				_rab[0] = _rb[0] - _ra[0]; // the vector connecting atoms alpha and beta
				_rab[1] = _rb[1] - _ra[1];
				_r = sqrt(pow(_rab[0], 2) + pow(_rab[1], 2)); // bond length

				// Get data for the bond connecting atoms alpha and beta
				_bond = (int)_BondList[_j];
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[2], _bond - 1,
					name_Potential));
				_z_p = prz[_bond - 1];

				// Compute dphi = \phi'
				_dphi = 0.5 * (_Potential[0] / _R) * (_r - _R - _z_p);

				// Compute force components
				_fa[0] = -_w*_dphi * _rab[0] / _r;
				_fa[1] = -_w*_dphi * _rab[1] / _r;

				// Allocate fa to Grad
#pragma omp atomic
				prGrad[2 * (_alpha - 1)] += _fa[0];
#pragma omp atomic
				prGrad[2 * (_alpha - 1) + 1] += _fa[1];

				// Allocate fa to Grad force acting on beta
#pragma omp atomic
				prGrad[2 * (_beta - 1)] -= _fa[0];
#pragma omp atomic
				prGrad[2 * (_beta - 1) + 1] -= _fa[1];
			}
		}
	}
}
