#include "qr.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

int main() {

	int N_rhs = 2;

	// QR decomposition example: algorithm (2) from arXiv:1710.09745
	// Given M, find Q and R such that QR = M, where:
	// M_i = fermion vector with block index i, i=0,..,N_rhs-1
	// Q_i = orthonormal fermion vector with block index i, i=0,..,N_rhs-1
	// R = N_rhsxN_rhs hermitian matrix

	// make dummy fermion vectors for this example 
	// with single complex double in place of each vector:
	_Complex double* Q;
	_Complex double* M;
	Q = (_Complex double*)malloc(N_rhs*sizeof(_Complex double));
	M = (_Complex double*)malloc(N_rhs*sizeof(_Complex double));
	for(int i=0; i<N_rhs; ++i) {
		M[i] = 0.345634 + 1.3245*i - 2.876*I*i; //some dummy values for M
	}

	// allocate space for hermitian matrices H and R
	// since they are hermitian we only need to store the lower triangular part
	double* H;
	double* R;
	int N_triangular = 0;
	for(int i=1; i<=N_rhs; ++i) {
		N_triangular += i;
	}
	H = (double*)malloc(N_triangular*2*sizeof(double));
	R = (double*)malloc(N_triangular*2*sizeof(double));

	// construct NxN hermitian matrix H_ij = (M_i, M_j)
	// for j<=i, from dot products of fermion vectors M
	int array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<=i; ++j) {
			_Complex double dot = conj(M[i]) * M[j];  // replace with vector dot product (M_i, M_j)
			H[array_index++] = creal(dot);
			H[array_index++] = cimag(dot);
		}
	}

	// call QR function defined in qr.h:
	// takes integer N_rhs and array of 2 N_rhs doubles H as input,
	// output goes in R: array of 2 N_rhs doubles
	qr(N_rhs, H, R);

	// Solve QR = M for Q by back substitution, where R is upper triangular
	array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		// Q_i = M_i
		Q[i] = M[i];
		for(int j=0; j<i; ++j) {
			// Q_i -= R(j,i) Q_j;
			double re = R[array_index++];
			double im = R[array_index++];
			Q[i] -= (re + I * im) * Q[j];
		}
		// Q_i /= R(i,i)
		Q[i] /= R[array_index++]; // R_ii is real
		array_index++; // skip complex part of R_ii
	}

	// Confirm that Q is orthonormal:
	printf("(Q,Q):\n");
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<N_rhs; ++j) {
			printf("%.8f\t", creal(conj(Q[i])*Q[j]));
		}
		printf("\n");
	}


	free(Q);
	free(M);
	free(H);
	free(R);
	return 0;
}