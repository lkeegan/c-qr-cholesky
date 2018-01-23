#include "qr.h"
#include <complex>
#include <iostream>
#include "Eigen/Dense"

extern "C" {

// Do QR decomposition of hermitian NxN matrix H, i.e. find Q and R such that H = QR
// where Q is orthonormal hermitian, and R is upper triangular
// Since H is hermitian, only the lower triangular elements are required which
// are stored as pairs of doubles representing the real and imaginary parts
// Moreover the diagonal elements are pure real so only require a single double
// Input format for H is assumed to be an array of doubles with row,col ordering:
// H_00_re, H_10_re, H_10_im, H_11_re, H_20_re, H_20_im, H_21_re, H_21_im, ...
void qr(int N, double* H_array, double* R_array) {
	// Construct lower triangular part of H from array of H elements
	Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(N, N);
	int array_index = 0;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<=i; ++j) {
			double re = H_array[array_index++];
			double im = 0.;
			if(i!=j) {
				// only off-diagonal elements are complex
				im = H_array[array_index++];
			}
			H(i, j) = std::complex<double>(re, im);
		}
	}

	// Find upper triangular R such that R^dag R = H = M^dag M
	// i.e. adjoint of cholesky decomposition L matrix: L L^dag = H
	Eigen::MatrixXcd R = H.llt().matrixL().adjoint();

	// Return TRANSPOSE of upper triangular R in same format as H:
	array_index = 0;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<=i; ++j) {
			R_array[array_index++] = R(j, i).real();
			if(i!=j) {
				// only off-diagonal elements are complex
				R_array[array_index++] = R(j, i).imag();
			}
		}
	}
}

// Find inverse of hermitian NxN matrix H
// Same format for H as for qr function:
// H_00_re, H_10_re, H_10_im, H_11_re, H_20_re, H_20_im, H_21_re, H_21_im, ...
void inverse(int N, double* H_array, double* inverseH_array) {
	// Construct lower triangular part of H from array of H elements
	Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(N, N);
	int array_index = 0;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<=i; ++j) {
			double re = H_array[array_index++];
			double im = 0.;
			if(i!=j) {
				// only off-diagonal elements are complex
				im = H_array[array_index++];
			}
			H(i, j) = std::complex<double>(re, im);
		}
	}

	// Find inverse of H via cholesky decomposition H = L L^dag
	// and solving H H^-1 = I for H^-1 back back/front substitution
	Eigen::MatrixXcd inverseH = H.llt().solve(Eigen::MatrixXcd::Identity(N, N));

	// Return H^-1 as array of doubles, same format as H:
	array_index = 0;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<=i; ++j) {
			inverseH_array[array_index++] = inverseH(i, j).real();
			if(i!=j) {
				// only off-diagonal elements are complex
				inverseH_array[array_index++] = inverseH(i, j).imag();
			}
		}
	}
}


} /* extern "C" */