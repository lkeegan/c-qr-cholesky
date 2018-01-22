#include "qr.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

// function to print lower triangular part of hermitian matrix
void print_lt_matrix(int N_rhs, double* A_lt) {
	int	array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<=i; ++j) {
			double re = A_lt[array_index++];
			double im = 0.;
			if(i!=j) {
				// only off-diagonal elements are complex
				im = A_lt[array_index++];
			}
			printf("(%.8f, %.8f)\t", re, im);
		}
		printf("\n");
	}
}

// function to print upper triangular matrix CURRENTLY WRONG ORDERING!
void print_ut_matrix(int N_rhs, double* A_ut) {
	int	array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<N_rhs; ++j) {
			double re = 0;
			double im = 0.;
			if (j==i) {
				re = A_ut[array_index++];
			} 
			else if (j>i) {
				re = A_ut[array_index++];
				im = A_ut[array_index++];
			}
			printf("(%.8f, %.8f)\t", re, im);
		}
		printf("\n");
	}
}

// function to multiply two hermitian matrices A and B,
// given only their lower triangular parts A_lt, B_lt,
// and print the resulting matrix
void print_lt_matrix_multiplication(int N_rhs, double* A_lt, double* B_lt) {
	// Allocate NxN empty complex matrices H, R, HR
	_Complex double** A = (_Complex double**)malloc(N_rhs*sizeof(_Complex double*));
	_Complex double** B = (_Complex double**)malloc(N_rhs*sizeof(_Complex double*));
	_Complex double** AB = (_Complex double**)malloc(N_rhs*sizeof(_Complex double*));
	for(int i=0; i<N_rhs; ++i) {
		A[i] = (_Complex double*)malloc(N_rhs*sizeof(_Complex double));
		B[i] = (_Complex double*)malloc(N_rhs*sizeof(_Complex double));
		AB[i] = (_Complex double*)malloc(N_rhs*sizeof(_Complex double));
		for(int j=0; j<N_rhs; ++j) {
			AB[i][j] = 0.0;
		}
	}
	// reconstruct A and B from lower triangular parts
	int array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<=i; ++j) {
			// real part:
			_Complex double A_ij = A_lt[array_index];
			_Complex double B_ij = B_lt[array_index++];
			if(i!=j) {
				// only off-diagonal elements are complex:
				A_ij += I*A_lt[array_index];
				B_ij += I*B_lt[array_index++];
			}
			A[i][j] = A_ij;
			A[j][i] = conj(A_ij);
			B[i][j] = B_ij;
			B[j][i] = conj(B_ij);
		}
	}
	// Do matrix multiplication
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<N_rhs; ++j) {
			for(int k=0; k<N_rhs; ++k) {
				AB[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	// Print result
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<N_rhs; ++j) {
			printf("(%.8e, %.8f)\t", creal(AB[i][j]), cimag(AB[i][j]));
		}
		printf("\n");
	}
	for(int i=0; i<N_rhs; ++i) {
		free(A[i]);
		free(B[i]);
		free(AB[i]);
	}
	free(A);
	free(B);
	free(AB);
}

// Inverse of hermitian matrix example
// Given NxN hermitian matrix H, find inverse R
// where RH = HR = I
void inverse_example(int N_rhs) {

	printf("\n\n\n=======================================\n");
	printf("  Inverse of hermitian matrix example\n");
	printf("=======================================\n\n");

	// allocate space for hermitian matrices H and R
	// require N_rhs doubles for diagonal real elements
	// plus 2 doubles for each complex element below the diagonal
	int N_triangular = N_rhs;
	for(int i=1; i<N_rhs; ++i) {
		N_triangular += 2*i;
	}
	double* H_lt = (double*)malloc(N_triangular*sizeof(double));
	double* R_lt = (double*)malloc(N_triangular*sizeof(double));

	// construct lower triangular part of NxN hermitian matrix H_ij
	// i.e. H_ij elements for j<=i
	int array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<=i; ++j) {
			_Complex double H_ij = 0.0043634*(i-j)*(1.1-I); // dummy values
			if(i==j) {
				H_ij += 1.234*(i+0.7); // dummy values
			}
			H_lt[array_index++] = creal(H_ij);
			if(i!=j) {
				// only off-diagonal elements are complex:
				H_lt[array_index++] = cimag(H_ij);
			}
		}
	}
	printf("Lower triangular part of hermitian matrix H:\n\n");
	print_lt_matrix(N_rhs, H_lt);

	// call inverse function defined in qr.h:
	// takes integer N_rhs and array of doubles containing 
	// lower triangular part of H as input, with ordering
	// H_00_re, H_00_im, H_10_re, H_10_im, H_11_re, H_11_im, H_20_re, H_20_im, ...
	// puts inverse of H into R in the same format.
	inverse(N_rhs, H_lt, R_lt);

	printf("\nLower triangular part of H^{-1}:\n\n");
	print_lt_matrix(N_rhs, R_lt);

	// Confirm that H*R = I
	printf("\nH * H^{-1} = Identity:\n\n");
	print_lt_matrix_multiplication(N_rhs, H_lt, R_lt);

	free(H_lt);
	free(R_lt);
}

_Complex double dot_product(int N_vol, _Complex double* A, _Complex double* B) {
	_Complex double dot = 0;
	for(int j=0; j<N_vol; ++j) {
		dot += conj(A[j])*B[j];
	}
	return dot;
}

// QR decomposition example: algorithm (2) from arXiv:1710.09745
// Given M, find Q and R such that QR = M, where:
// M_i = fermion vector with N_vol complex components, block index i: i=0,..,N_rhs-1
// Q_i = orthonormal version of M
// R = N_rhs x N_rhs hermitian matrix
void qr_example(int N_rhs, int N_vol) {
	printf("\n\n\n=======================================\n");
	printf("       QR decomposition example\n");
	printf("=======================================\n\n");

	// make dummy fermion vectors M for this example 
	_Complex double** M = (_Complex double**)malloc(N_rhs*sizeof(_Complex double*));
	for(int i=0; i<N_rhs; ++i) {
		M[i] = (_Complex double*)malloc(N_vol*sizeof(_Complex double));
		for(int j=0; j<N_vol; ++j) {
			M[i][j] = 0.02*i + 0.08*j + 0.70*I*i - 0.058*I*j;
			if(j==i) {
				M[i][j] += 1.0;
			}
		}
	}
	// allocate fermion vectors Q
	_Complex double** Q = (_Complex double**)malloc(N_rhs*sizeof(_Complex double*));
	for(int i=0; i<N_rhs; ++i) {
		Q[i] = (_Complex double*)malloc(N_vol*sizeof(_Complex double));
	}

	// allocate space for hermitian matrix H and upper triangular R
	// require N_rhs doubles for diagonal real elements
	// plus 2 doubles for each complex element above the diagonal
	int N_triangular = N_rhs;
	for(int i=1; i<N_rhs; ++i) {
		N_triangular += 2*i;
	}
	double* H_lt = (double*)malloc(N_triangular*sizeof(double));
	double* R_ut = (double*)malloc(N_triangular*sizeof(double));

	// construct NxN hermitian matrix H_ij = (M_i, M_j)
	// for j<=i, from dot products of fermion vectors M
	int array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<=i; ++j) {
			// vector dot product H_ij = (M_i, M_j)
			_Complex double H_ij = dot_product(N_vol, M[i], M[j]);
			H_lt[array_index++] = creal(H_ij);
			if(i!=j) {
				// only off-diagonal elements are complex:
				H_lt[array_index++] = cimag(H_ij);
			}
		}
	}
	printf("Lower triangular part of hermitian matrix H:\n\n");
	print_lt_matrix(N_rhs, H_lt);

	// call QR function defined in qr.h:
	// takes integer N_rhs and array of doubles containing 
	// lower triangular part of H as input, with ordering
	// H_00_re, H_00_im, H_10_re, H_10_im, H_11_re, H_11_im, H_20_re, H_20_im, ...
	// puts **TODO** output into R in the same format, where QR = H
	qr(N_rhs, H_lt, R_ut);
	printf("Upper triangular matrix R:\n\n");
	print_ut_matrix(N_rhs, R_ut);

	// Solve QR = M for Q by back substitution, where R is upper triangular
	array_index = 0;
	for(int i=0; i<N_rhs; ++i) {
		// Q_i = M_i
		for(int k=0; k<N_vol; ++k) {
			Q[i][k] = M[i][k];
		}
		for(int j=0; j<i; ++j) {
			// Q_i -= R(j,i) Q_j;
			double re = R_ut[array_index++];
			double im = R_ut[array_index++];
			for(int k=0; k<N_vol; ++k) {
				Q[i][k] -= (re + I * im) * Q[j][k];
			}
		}
		// Q_i /= R(i,i)
		double re = R_ut[array_index++]; // R_ii is real
		for(int k=0; k<N_vol; ++k) {
			Q[i][k] /= re;
		}
	}

	// Confirm that Q vectors are orthonormal:
	printf("(Q,Q):\n");
	for(int i=0; i<N_rhs; ++i) {
		for(int j=0; j<N_rhs; ++j) {
			_Complex double dot = dot_product(N_vol, Q[i], Q[j]);
			printf("(%.8f, %.8f)\t", creal(dot), cimag(dot));
		}
		printf("\n");
	}

	free(H_lt);
	free(R_ut);
	for(int i=0; i<N_rhs; ++i) {
		free(Q[i]);
		free(M[i]);
	}
	free(Q);
	free(M);
}

int main() {

	// example inversion of 6x6 hermitian matrix:
	inverse_example(6);

	// example thinQR for n_RHS = 5, n_VOL = 400
	qr_example(5, 400);

	return 0;
}