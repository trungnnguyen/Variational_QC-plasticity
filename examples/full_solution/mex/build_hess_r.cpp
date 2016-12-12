// build_hess_r: assembly the Hessian of internal energy w.r.t. r

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
	int NBonds = (int)mxGetNumberOfElements(prhs[1]);

	// Allocate outputs
	nlhs = 3;
	plhs[0] = mxCreateDoubleMatrix(16, NBonds, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(16, NBonds, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(16, NBonds, mxREAL);
	double* prI = mxGetPr(plhs[0]);
	double* prJ = mxGetPr(plhs[1]);
	double* prS = mxGetPr(plhs[2]);

	// Compute the Hessian and populate prI, prJ, prS
	int i;
	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[1], "Potential");

#pragma omp parallel
	{
		int _j, _m, _n, _alpha, _beta, _bond;
		double _R, _r, _dphi, _ddphi, _z_p;
		double *_Ra, *_Rb, *_NeighbourList, *_BondList, *_Potential;
		double _ra[2], _rb[2], _rab[2];
		double _Kbond[4][4]; // local stiffness matrix for a single bond
		int _Ln[4]; // vector of localization numbers
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
				_z_p = prz[_bond - 1];

				// Compute dphi = \phi', and ddphi = \phi''
				_dphi = 0.5*(_Potential[0] / _R)*(_r - _R - _z_p);
				_ddphi = 0.5*_Potential[0] / _R;

				// Construct Kbond - hessian matrix of a single bond
				// gamma = alpha, delta = alpha
				_Kbond[0][0] = (_dphi / _r * 1 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[0] * _rab[0]));
				_Kbond[0][1] = (_dphi / _r * 0 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[0] * _rab[1]));
				_Kbond[1][0] = (_dphi / _r * 0 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[1] * _rab[0]));
				_Kbond[1][1] = (_dphi / _r * 1 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[1] * _rab[1]));
				// gamma = beta, delta = beta
				_Kbond[2][2] = _Kbond[0][0];
				_Kbond[2][3] = _Kbond[0][1];
				_Kbond[3][2] = _Kbond[1][0];
				_Kbond[3][3] = _Kbond[1][1];
				// gamma = alpha, delta = beta
				_Kbond[0][2] = -(_dphi / _r * 1 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[0] * _rab[0]));
				_Kbond[0][3] = -(_dphi / _r * 0 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[0] * _rab[1]));
				_Kbond[1][2] = -(_dphi / _r * 0 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[1] * _rab[0]));
				_Kbond[1][3] = -(_dphi / _r * 1 + (_ddphi / (_r*_r) -
					_dphi / (_r*_r*_r))*(_rab[1] * _rab[1]));
				// gamma = beta, delta = alpha
				_Kbond[2][0] = _Kbond[0][2];
				_Kbond[2][1] = _Kbond[0][3];
				_Kbond[3][0] = _Kbond[1][2];
				_Kbond[3][1] = _Kbond[1][3];

				// Allocate Kbond into prI, prJ, prS
				_Ln[0] = 2 * _alpha - 1; // indexing for MATLAB
				_Ln[1] = 2 * _alpha;
				_Ln[2] = 2 * _beta - 1;
				_Ln[3] = 2 * _beta;

#pragma omp critical(ALLOCATE)
				{
					for (_n = 0; _n < 4; _n++){
						for (_m = 0; _m < 4; _m++){
							prI[(_bond - 1) * 16 + _n * 4 + _m] = _Ln[_m]; // note that bonds are counted twice
							prJ[(_bond - 1) * 16 + _n * 4 + _m] = _Ln[_n];
							prS[(_bond - 1) * 16 + _n * 4 + _m] +=
								_Kbond[_m][_n];
						}
					}
				}
			}
		}
	}
}
