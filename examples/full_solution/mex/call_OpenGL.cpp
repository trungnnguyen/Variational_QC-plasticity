// call_OpenGL: creates shared memory arrays for all the necessary data and
// calls draw_OpenGL.exe postprocessing tool

#include <stdio.h>
#include <windows.h>
#include "mex.h"

/************************ Main program ************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[]){

	// Test for eleven input arguments
	if (nrhs != 11){
		mexErrMsgTxt("Eleven input arguments required.");
	}

	// Create data to share
	// Atoms
	// Satoms
	int i, j;
	double *prID;
	int NsAtoms = (int)mxGetNumberOfElements(prhs[1]);
	INT* satoms = new INT[NsAtoms + 2];
	int name_ID = mxGetFieldNumber(prhs[1], "ID");
	satoms[0] = NsAtoms; // number of rows
	satoms[1] = 1; // number of cols
	for (i = 0; i < satoms[0]; i++){
		prID = mxGetPr(mxGetFieldByNumber(prhs[1], i, name_ID));
		for (j = 0; j < satoms[1]; j++){
			satoms[(i + NsAtoms*j) + 2] = (INT)prID[0];
		}
	}

	// AtomsQC
	double *prAtomsQC = mxGetPr(prhs[2]);
	int NrepAtoms = (int)mxGetM(prhs[2]);
	INT* atomsQC = new INT[NrepAtoms + 2];
	atomsQC[0] = NrepAtoms; // number of rows
	atomsQC[1] = 1; // number of cols
	for (i = 0; i < atomsQC[0]; i++){
		for (j = 0; j < atomsQC[1]; j++){
			atomsQC[(i + NrepAtoms*j) + 2] =
				(INT)prAtomsQC[i + NrepAtoms*j];
		}
	}

	// Bonds - columnwise: Atoms alpha, beta, Potential
	double *prAtoms, *prPotential;
	int name_Atoms = mxGetFieldNumber(prhs[3], "Atoms");
	int name_Potential = mxGetFieldNumber(prhs[3], "Potential");
	int NBonds = (int)mxGetNumberOfElements(prhs[3]);
	FLOAT* bonds = new FLOAT[NBonds * 8 + 2];
	bonds[0] = (FLOAT)NBonds; // number of rows
	bonds[1] = (FLOAT)(2 + 6); // number of cols
	for (i = 0; i < bonds[0]; i++){
		prAtoms = mxGetPr(mxGetFieldByNumber(prhs[3], i, name_Atoms));
		prPotential = mxGetPr(mxGetFieldByNumber(prhs[3], i,
			name_Potential));
		for (j = 0; j < 2; j++){
			bonds[(i + NBonds*j) + 2] = (FLOAT)prAtoms[j];
		}
		for (j = 2; j < 8; j++){
			bonds[(i + NBonds*j) + 2] = (FLOAT)prPotential[j - 2];
		}
	}

	// BondsQC
	double *prBondsQC = mxGetPr(prhs[4]);
	int NBondsQC = (int)mxGetM(prhs[4]);
	INT* bondsQC = new INT[NBondsQC + 2];
	bondsQC[0] = NBondsQC; // number of rows
	bondsQC[1] = 1; // number of cols
	for (i = 0; i < bondsQC[0]; i++){
		for (j = 0; j < bondsQC[1]; j++){
			bondsQC[(i + NBondsQC*j) + 2] = (INT)prBondsQC[i + NBondsQC*j];
		}
	}

	// Triangles - row-wise stores VertexAtoms
	double *prVertexAtoms;
	int Ntriangles = (int)mxGetNumberOfElements(prhs[5]);
	INT* triangles = new INT[Ntriangles * 3 + 2];
	int name_VertexAtoms = mxGetFieldNumber(prhs[5], "VertexAtoms");
	triangles[0] = Ntriangles; // number of rows
	triangles[1] = 3; // number of cols
	for (i = 0; i < triangles[0]; i++){
		prVertexAtoms = mxGetPr(mxGetFieldByNumber(prhs[5], i,
			name_VertexAtoms));
		for (j = 0; j < triangles[1]; j++){
			triangles[(i + Ntriangles*j) + 2] = (INT)prVertexAtoms[j];
		}
	}

	// R
	double *prR = mxGetPr(prhs[6]);
	int MR = (int)mxGetM(prhs[6]);
	int NR = (int)mxGetN(prhs[6]);
	FLOAT* R = new FLOAT[MR*NR + 2];
	R[0] = (FLOAT)MR; // number of rows
	R[1] = (FLOAT)NR; // number of cols
	for (i = 0; i < R[0]; i++){
		for (j = 0; j < R[1]; j++){
			R[(i + MR*j) + 2] = (FLOAT)prR[i + MR*j];
		}
	}

	// Z
	double *prZ = mxGetPr(prhs[7]);
	int MZ = (int)mxGetM(prhs[7]);
	int NZ = (int)mxGetN(prhs[7]);
	FLOAT* Z = new FLOAT[MZ*NZ + 2];
	Z[0] = (FLOAT)MZ; // number of rows
	Z[1] = (FLOAT)NZ; // number of cols
	for (i = 0; i < Z[0]; i++){
		for (j = 0; j < Z[1]; j++){
			Z[(i + MZ*j) + 2] = (FLOAT)prZ[i + MZ*j];
		}
	}

	// K
	double *prK = mxGetPr(prhs[8]);
	int MK = (int)mxGetM(prhs[8]);
	int NK = (int)mxGetN(prhs[8]);
	FLOAT* K = new FLOAT[MK*NK + 2];
	K[0] = (FLOAT)MK; // number of rows
	K[1] = (FLOAT)NK; // number of cols
	for (i = 0; i < K[0]; i++){
		for (j = 0; j < K[1]; j++){
			K[(i + MK*j) + 2] = (FLOAT)prK[i + MK*j];
		}
	}

	// Time
	double *prTime = mxGetPr(prhs[9]);
	int MTime = (int)mxGetM(prhs[9]);
	int NTime = (int)mxGetN(prhs[9]);
	FLOAT* Time = new FLOAT[MTime*NTime + 2];
	Time[0] = (FLOAT)MTime; // number of rows
	Time[1] = (FLOAT)NTime; // number of cols
	for (i = 0; i < Time[0]; i++){
		for (j = 0; j < Time[1]; j++){
			Time[(i + MTime*j) + 2] = (FLOAT)prTime[i + MTime*j];
		}
	}

	// Size
	double *prSize = mxGetPr(prhs[10]);
	int MSize = (int)mxGetM(prhs[10]);
	int NSize = (int)mxGetN(prhs[10]);
	FLOAT* Size = new FLOAT[MSize*NSize + 2];
	Size[0] = (FLOAT)MSize; // number of rows
	Size[1] = (FLOAT)NSize; // number of cols
	for (i = 0; i < Size[0]; i++){
		for (j = 0; j < Size[1]; j++){
			Size[(i + MSize*j) + 2] = (FLOAT)prSize[i + MSize*j];
		}
	}

	// Create named shared memory
	HANDLE hsatoms, hatomsQC, hbonds, hbondsQC, htriangles, hR, hZ, hK,
		hTime, hSize;
	INT *psatoms, *patomsQC, *ptriangles, *pbondsQC;
	FLOAT *pbonds, *pR, *pZ, *pK, *pTime, *pSize;
	hsatoms = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
		PAGE_READWRITE, 0, (NsAtoms + 2)*sizeof(INT), "satoms");
	hatomsQC = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
		PAGE_READWRITE, 0, (NrepAtoms + 2)*sizeof(INT), "atomsQC");
	hbonds = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
		PAGE_READWRITE, 0, (NBonds * 8 + 2)*sizeof(FLOAT), "bonds");
	hbondsQC = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
		PAGE_READWRITE, 0, (NBondsQC + 2)*sizeof(INT), "bondsQC");
	htriangles = CreateFileMapping(INVALID_HANDLE_VALUE, NULL,
		PAGE_READWRITE, 0, (Ntriangles * 3 + 2)*sizeof(INT), "triangles");
	hR = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
		0, (MR*NR + 2)*sizeof(FLOAT), "R");
	hZ = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
		0, (MZ*NZ + 2)*sizeof(FLOAT), "Z");
	hK = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
		0, (MK*NK + 2)*sizeof(FLOAT), "K");
	hTime = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
		0, (MTime*NTime + 2)*sizeof(FLOAT), "Time");
	hSize = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE,
		0, (MSize*NSize + 2)*sizeof(FLOAT), "size");
	if ((hsatoms == NULL) || (hatomsQC == NULL) || (hbonds == NULL) ||
		(hbondsQC == NULL) || (htriangles == NULL) || (hR == NULL) ||
		(hZ == NULL) || (hK == NULL) || (hTime == NULL) ||
		(hSize == NULL)){
		delete[] satoms;
		delete[] atomsQC;
		delete[] bonds;
		delete[] bondsQC;
		delete[] triangles;
		delete[] R;
		delete[] Z;
		delete[] K;
		delete[] Time;
		delete[] Size;

		mexErrMsgTxt("Could not create file mapping objects.");
	}
	psatoms = (INT*)MapViewOfFile(hsatoms, FILE_MAP_ALL_ACCESS, 0, 0,
		(NsAtoms + 2)*sizeof(INT));
	patomsQC = (INT*)MapViewOfFile(hatomsQC, FILE_MAP_ALL_ACCESS, 0, 0,
		(NrepAtoms + 2)*sizeof(INT));
	pbonds = (FLOAT*)MapViewOfFile(hbonds, FILE_MAP_ALL_ACCESS, 0, 0,
		(NBonds * 8 + 2)*sizeof(FLOAT));
	pbondsQC = (INT*)MapViewOfFile(hbondsQC, FILE_MAP_ALL_ACCESS, 0, 0,
		(NBondsQC + 2)*sizeof(INT));
	ptriangles = (INT*)MapViewOfFile(htriangles, FILE_MAP_ALL_ACCESS, 0, 0,
		(Ntriangles * 3 + 2)*sizeof(INT));
	pR = (FLOAT*)MapViewOfFile(hR, FILE_MAP_ALL_ACCESS, 0, 0,
		(MR*NR + 2)*sizeof(FLOAT));
	pZ = (FLOAT*)MapViewOfFile(hZ, FILE_MAP_ALL_ACCESS, 0, 0,
		(MZ*NZ + 2)*sizeof(FLOAT));
	pK = (FLOAT*)MapViewOfFile(hK, FILE_MAP_ALL_ACCESS, 0, 0,
		(MK*NK + 2)*sizeof(FLOAT));
	pTime = (FLOAT*)MapViewOfFile(hTime, FILE_MAP_ALL_ACCESS, 0, 0,
		(MTime*NTime + 2)*sizeof(FLOAT));
	pSize = (FLOAT*)MapViewOfFile(hSize, FILE_MAP_ALL_ACCESS, 0, 0,
		(MSize*NSize + 2)*sizeof(FLOAT));
	if ((psatoms == NULL) || (patomsQC == NULL) || (pbonds == NULL) ||
		(pbondsQC == NULL) || (ptriangles == NULL) || (pR == NULL) ||
		(pZ == NULL) || (pK == NULL) || (pTime == NULL) ||
		(pSize == NULL)){
		CloseHandle(hsatoms);
		CloseHandle(hatomsQC);
		CloseHandle(hbonds);
		CloseHandle(hbondsQC);
		CloseHandle(htriangles);
		CloseHandle(hR);
		CloseHandle(hZ);
		CloseHandle(hK);
		CloseHandle(hTime);
		CloseHandle(hSize);

		delete[] satoms;
		delete[] atomsQC;
		delete[] bonds;
		delete[] bondsQC;
		delete[] triangles;
		delete[] R;
		delete[] Z;
		delete[] K;
		delete[] Time;
		delete[] Size;

		mexErrMsgTxt("Could not map view of files.");
	}
	CopyMemory((PVOID)psatoms, satoms, (NsAtoms + 2)*sizeof(INT));
	CopyMemory((PVOID)patomsQC, atomsQC, (NrepAtoms + 2)*sizeof(INT));
	CopyMemory((PVOID)pbonds, bonds, (NBonds * 8 + 2)*sizeof(FLOAT));
	CopyMemory((PVOID)pbondsQC, bondsQC, (NBondsQC + 2)*sizeof(INT));
	CopyMemory((PVOID)ptriangles, triangles,
		(Ntriangles * 3 + 2)*sizeof(INT));
	CopyMemory((PVOID)pR, R, (MR*NR + 2)*sizeof(FLOAT));
	CopyMemory((PVOID)pZ, Z, (MZ*NZ + 2)*sizeof(FLOAT));
	CopyMemory((PVOID)pK, K, (MK*NK + 2)*sizeof(FLOAT));
	CopyMemory((PVOID)pTime, Time, (MTime*NTime + 2)*sizeof(FLOAT));
	CopyMemory((PVOID)pSize, Size, (MSize*NSize + 2)*sizeof(FLOAT));

	// Call the second process - draw_OpenGL.exe
	printf("\nOpenGL called.\n");
	system("draw_OpenGL.exe");

	// Unmap buffers
	UnmapViewOfFile(psatoms);
	UnmapViewOfFile(patomsQC);
	UnmapViewOfFile(pbonds);
	UnmapViewOfFile(pbondsQC);
	UnmapViewOfFile(ptriangles);
	UnmapViewOfFile(pR);
	UnmapViewOfFile(pZ);
	UnmapViewOfFile(pK);
	UnmapViewOfFile(pTime);
	UnmapViewOfFile(pSize);
	CloseHandle(hsatoms);
	CloseHandle(hatomsQC);
	CloseHandle(hbonds);
	CloseHandle(hbondsQC);
	CloseHandle(htriangles);
	CloseHandle(hR);
	CloseHandle(hZ);
	CloseHandle(hK);
	CloseHandle(hTime);
	CloseHandle(hSize);

	delete[] satoms;
	delete[] atomsQC;
	delete[] bonds;
	delete[] bondsQC;
	delete[] triangles;
	delete[] R;
	delete[] Z;
	delete[] K;
	delete[] Time;
	delete[] Size;

	printf("OpenGL closed.\n");
}
