// return_mapping_QC: compute Z using return-mapping algorithm for all 
// sampling bonds

#include <stdio.h>
#include <math.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mex.h"
#include "matrix.h"

/******************* Definition of constants ******************/
const int MAXITER = 100; // maximum No. of Newton iterations

/******************** Function declarations *******************/
int sgn(double x); // signum function sgn(x)

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for seven input arguments
	if (nrhs != 7){
		mexErrMsgTxt("Seven input arguments required.");
	}

	// Get the input data
	double* prsampbondsID = mxGetPr(prhs[2]);
	double* prr = mxGetPr(prhs[3]);
	double* prz = mxGetPr(prhs[4]);
	double* prk = mxGetPr(prhs[5]);
	double* prTOL_z = mxGetPr(prhs[6]);
	double TOL_z = prTOL_z[0];

	// Minimize bond-wise energies
	nlhs = 1;
	int i;
	int NBonds = (int)mxGetNumberOfElements(prhs[1]);
	int NsBonds = (int)mxGetNumberOfElements(prhs[2]);
	plhs[0] = mxCreateDoubleMatrix(NBonds, 1, mxREAL);
	double* z_p = mxGetPr(plhs[0]); // current plastic elongations
	int name_R = mxGetFieldNumber(prhs[0], "R");
	int name_Atoms = mxGetFieldNumber(prhs[1], "Atoms");
	int name_Potential = mxGetFieldNumber(prhs[1], "Potential");

#pragma omp parallel
	{
		int _alpha, _beta, _bond;
		double _R, _r, _Zp, _Zc, _ftri, _Yieldtri, _tempz_p, _dZp;
		double *_Ra, *_Rb, *_Potential, *_Atoms, _ra[2], _rb[2];

#pragma omp for
		for (i = 0; i < NsBonds; i++){ // loop over all sampling bonds
			_bond = (int)prsampbondsID[i]; // bond ID, MATLAB indexing

			// Get data for atoms alpha and beta connected by the bond
			_Atoms = mxGetPr(mxGetFieldByNumber(prhs[1], _bond - 1, name_Atoms));
			_alpha = (int)_Atoms[0]; // atom alpha, MATLAB indexing
			_beta = (int)_Atoms[1]; // atom beta, MATLAB indexing
			_Ra = mxGetPr(mxGetFieldByNumber(prhs[0], _alpha - 1, name_R));
			_Rb = mxGetPr(mxGetFieldByNumber(prhs[0], _beta - 1, name_R));
			_ra[0] = prr[2 * (_alpha - 1)];
			_ra[1] = prr[2 * (_alpha - 1) + 1];
			_rb[0] = prr[2 * (_beta - 1)];
			_rb[1] = prr[2 * (_beta - 1) + 1];
			_R = sqrt(pow(_Rb[0] - _Ra[0], 2) + pow(_Rb[1] - _Ra[1], 2));
			_r = sqrt(pow(_rb[0] - _ra[0], 2) + pow(_rb[1] - _ra[1], 2)); // bond length

			// Get data for the bond between atoms alpha and beta
			_Potential = mxGetPr(mxGetFieldByNumber(prhs[1], _bond - 1,
				name_Potential));
			_Zp = prz[_bond - 1]; // previous time-step plastic elongation
			_Zc = prk[_bond - 1]; // previous time-step cumulative plastic elongation

			// Nonlinear hardeing
			_ftri = (_Potential[0] / _R)*(_r - _R - _Zp); // compute elastic trial force
			_Yieldtri = fabs(_ftri) -
				_Potential[2] * (1 + _Potential[1] *
				pow(_Zc / _R, _Potential[3]));

			// Elastic step
			if (_Yieldtri <= 0){
				_tempz_p = _Zp;
			}

			// Plastic step - return-mapping
			else{
				double _eps = TOL_z + 1;
				int _Niter = 0;
				double _x = 1e-14; // perturb searched value
				double _dx, _gx, _dgx, _normEps;
				while (_eps > TOL_z){ // 1-D Newton iteration
					_Niter++;
					_gx = fabs(_ftri) - _x*(_Potential[0] / _R) -
						_Potential[2] * (1 + _Potential[1] *
						pow((_Zc + _x) / _R, _Potential[3]));
					_dgx = -(_Potential[0] / _R) - _Potential[2] *
						_Potential[1] * pow((_Zc + _x) / _R,
						_Potential[3] - 1) / _R;
					_dx = -_gx / _dgx;
					_dx = std::max(double(0), _dx); // allow only for positive increments
					_x += _dx;
					if (_Niter == 1){
						_normEps = fabs(_dx);
					}
					_eps = fabs(_dx) / _normEps;
					if (_dx == 0){
						_eps = 0;
					}
					if (_Niter > MAXITER){
						printf("%u ReturnMapping iterations exceeded.\n", MAXITER);
						break;
					}
				}
				_dZp = _x;
				_tempz_p = _Zp + _dZp*sgn(_ftri);
			}

#pragma omp atomic
			z_p[_bond - 1] += _tempz_p;
		}
	}
}

/******************** Function definitions ********************/
// Signum function
int sgn(double x){
	return (x > 0) - (x < 0);
}
