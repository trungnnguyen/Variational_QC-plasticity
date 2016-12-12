// build_grad_r: assembly the gradient of the internal energy w.r.t. r

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

	// Test for four input arguments
	if (nrhs != 4){
		mexErrMsgTxt("Four input arguments required.");
	}

	// Get the input data
	double* prr = mxGetPr(prhs[2]);
	double* prz = mxGetPr(prhs[3]);
	int NAtoms = (int)mxGetNumberOfElements(prhs[0]);

	// Compute the gradient and send out the data back to MATLAB
	nlhs = 1;
	mxArray* Grad = mxCreateDoubleMatrix(2 * NAtoms, 1, mxREAL);
	double* prGrad = mxGetPr(Grad);
	plhs[0] = Grad;
	int i;

	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[1], "Potential");

#pragma omp parallel
	{
		int _j, _alpha, _beta, _bond;
		double _R, _r, _dphi, _z_p;
		double *_Ra, *_Rb, *_NeighbourList, *_BondList, *_Potential;
		double _ra[2], _rb[2], _rab[2], _fa[2];
		mxArray *_pNeighbourList; // pointer to NeighbourList array

#pragma omp for
		for (i = 0; i < NAtoms; i++){ // loop over all atoms
			_alpha = i + 1; // atom alpha ID, MATLAB indexing

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

				// Get data for atom beta
				_Rb = mxGetPr(mxGetFieldByNumber(prhs[0], _beta - 1,
					name_R));
				_rb[0] = prr[2 * (_beta - 1)];
				_rb[1] = prr[2 * (_beta - 1) + 1];
				_R = sqrt(pow(_Rb[0] - _Ra[0], 2) +
					pow(_Rb[1] - _Ra[1], 2)); // the initial bond length
				_rab[0] = _rb[0] - _ra[0]; // vector connecting atoms alpha and beta
				_rab[1] = _rb[1] - _ra[1];
				_r = sqrt(pow(_rab[0], 2) + pow(_rab[1], 2)); // bond length

				// Get data for the bond connecting atoms alpha and beta
				_bond = (int)_BondList[_j]; // MATLAB indexing
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[1], _bond - 1,
					name_Potential));
				_z_p = prz[_bond - 1]; // plastic deformation z_p

				// Compute dphi = \phi'
				_dphi = 0.5 * (_Potential[0] / _R) * (_r - _R - _z_p);

				// Compute force components
				_fa[0] = -_dphi * _rab[0] / _r;
				_fa[1] = -_dphi * _rab[1] / _r;

				// Allocate fa to Grad
#pragma omp atomic
				prGrad[2 * (_alpha - 1)] += _fa[0];
#pragma omp atomic
				prGrad[2 * (_alpha - 1) + 1] += _fa[1];

				// Allocate to Grad the force acting on atom beta
#pragma omp atomic
				prGrad[2 * (_beta - 1)] -= _fa[0];
#pragma omp atomic
				prGrad[2 * (_beta - 1) + 1] -= _fa[1];
			}
		}
	}
}
