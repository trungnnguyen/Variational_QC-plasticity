// build_hess_r_QC: assembly the Hessian of approximate internal 
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
	int NsAtoms = (int)mxGetNumberOfElements(prhs[1]);

	// Allocate outputs
	nlhs = 3;
	plhs[0] = mxCreateDoubleMatrix(16, 8 * NsAtoms, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(16, 8 * NsAtoms, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(16, 8 * NsAtoms, mxREAL);
	double* prI = mxGetPr(plhs[0]);
	double* prJ = mxGetPr(plhs[1]);
	double* prS = mxGetPr(plhs[2]);

	// Compute the Hessian and populate prI, prJ, prS
	int i;
	int counter = 0;
	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_NeighbourList = mxGetFieldNumber(prhs[0], "NeighbourList");
	int name_BondList = mxGetFieldNumber(prhs[0], "BondList");
	int name_Potential = mxGetFieldNumber(prhs[2], "Potential");
	int name_ID = mxGetFieldNumber(prhs[1], "ID");
	int name_w = mxGetFieldNumber(prhs[1], "w");

#pragma omp parallel
	{
		int _j, _m, _n, _alpha, _beta, _bond;
		double *_Ra, *_Rb, *_NeighbourList, *_BondList, *_prID,
			*_Potential, *_prw;
		double _ra[2], _rb[2], _rab[2], _R, _r, _dphi, _ddphi, _z_p, _w;
		double _Kbond[4][4]; // local stiffness matrix for a single bond
		int _Ln[4]; // vector of code numbers
		int* _prI = new int[16 * 8 * NsAtoms]; // thread-private row indices
		int* _prJ = new int[16 * 8 * NsAtoms]; // thread-private column indices
		double* _prS = new double[16 * 8 * NsAtoms]; // thread-private data
		int _counter = 0; // thread-private bond counter
		mxArray *_pNeighbourList;

#pragma omp for
		for (i = 0; i < NsAtoms; i++){ // loop over all sampling atoms
			_prw = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_w));
			_prID = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_ID));
			_w = _prw[0];
			_alpha = (int)_prID[0]; // atom alpha ID, MATLAB indexing

			// Get data for a sampling atom alpha
			_Ra = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1, name_R));
			_pNeighbourList = mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_NeighbourList);
			_NeighbourList = mxGetPr(_pNeighbourList);
			_BondList = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1,
				name_BondList));
			_ra[0] = prr[2 * (_alpha - 1)];
			_ra[1] = prr[2 * (_alpha - 1) + 1];
			for (_j = 0; _j < mxGetN(_pNeighbourList); _j++){ // loop over all the nearest neighbours
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
				_Potential = mxGetPr(mxGetFieldByNumber(prhs[2], _bond - 1,
					name_Potential));
				_z_p = prz[_bond - 1];

				// Compute dphi = \phi', and ddphi = \phi''
				_dphi = 0.5*(_Potential[0] / _R)*(_r - _R - _z_p);
				_ddphi = 0.5*_Potential[0] / _R;

				// Construct _Kbond - hessian matrix of a single bond
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

				// Allocate _Kbond into _prI, _prJ, and _prS 
				_Ln[0] = 2 * _alpha - 1; // indexing for MATLAB
				_Ln[1] = 2 * _alpha;
				_Ln[2] = 2 * _beta - 1;
				_Ln[3] = 2 * _beta;
				for (_n = 0; _n < 4; _n++){
					for (_m = 0; _m < 4; _m++){
						_prI[_counter * 16 + _n * 4 + _m] = _Ln[_m];
						_prJ[_counter * 16 + _n * 4 + _m] = _Ln[_n];
						_prS[_counter * 16 + _n * 4 + _m] =
							_w*_Kbond[_m][_n]; // multiply by the weight factor
					}
				}
				_counter++;
			}
		}

		// Collect _prI, _prJ, _prS from all threads
#pragma omp critical(ALLOCATE)
		{
			for (_j = 0; _j < _counter; _j++){
				for (_m = 0; _m < 16; _m++){
					prI[(counter + _j) * 16 + _m] = _prI[_j * 16 + _m];
					prJ[(counter + _j) * 16 + _m] = _prJ[_j * 16 + _m];
					prS[(counter + _j) * 16 + _m] = _prS[_j * 16 + _m];
				}
			}
			counter += _counter;
			delete[] _prI;
			delete[] _prJ;
			delete[] _prS;
		}
	}

	// Resize the outputs
	mxSetN(plhs[0], counter);
	mxSetN(plhs[1], counter);
	mxSetN(plhs[2], counter);
}
