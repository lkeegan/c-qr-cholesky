#include "qr.h"
#include <complex>
#include <iostream>
#include "Eigen/Dense"

extern "C" {

// Do QR decomposition of hermitian NxN matrix H, i.e. H = QR, output R
// Since H and R are hermitian, only the lower triangular elements are required
// Input format for H and R is assumed to be an array of doubles:
// H_00_re, H_00_im, H_10_re, H_10_im, H_11_re, H_11_im, H_20_re, H_20_im, ...
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
	std::cout << "\nSupplied H matrix:\n\n" << H << std::endl;

	// Find upper triangular R such that R^dag R = H = M^dag M
	// i.e. adjoint of cholesky decomposition L matrix: L L^dag = H
	Eigen::MatrixXcd R = H.llt().matrixL().adjoint();
	std::cout << "\nCalculated R matrix:\n\n" << R << std::endl;

	std::cout << H.llt().solve(Eigen::MatrixXcd::Identity(N, N)) - H.inverse() << std::endl;
	// Return upper triangular R as array of real doubles:
	array_index = 0;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<=i; ++j) {
			std::cout << j << "," << i << R(j,i) << std::endl;
			R_array[array_index++] = R(j, i).real();
			if(i!=j) {
				// only off-diagonal elements are complex
				R_array[array_index++] = R(j, i).imag();
			}
		}
	}
}

// Find inverse R of hermitian NxN matrix H, i.e. HR = RH = I, output R
// Since H is hermitian, only the lower triangular elements are required
// Format for  H and R is an array of doubles:
// H_00_re, H_00_im, H_10_re, H_10_im, H_11_re, H_11_im, H_20_re, H_20_im, ...
void inverse(int N, double* H_array, double* R_array) {
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
	Eigen::MatrixXcd R = H.llt().solve(Eigen::MatrixXcd::Identity(N, N));

	// Return R as array of doubles:
	array_index = 0;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<=i; ++j) {
			R_array[array_index++] = R(i, j).real();
			if(i!=j) {
				// only off-diagonal elements are complex
				R_array[array_index++] = R(i, j).imag();
			}
		}
	}
}


} /* extern "C" */