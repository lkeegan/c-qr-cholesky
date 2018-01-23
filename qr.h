#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// Do QR decomposition of hermitian NxN matrix H, i.e. find Q and R such that H = QR
// where Q is orthonormal hermitian, and R is upper triangular
// Since H is hermitian, only the lower triangular elements are required which
// are stored as pairs of doubles representing the real and imaginary parts
// Moreover the diagonal elements are pure real so only require a single double
// Input format for H is assumed to be an array of doubles with row,col ordering:
// H_00_re, H_10_re, H_10_im, H_11_re, H_20_re, H_20_im, H_21_re, H_21_im, ...
// Returns transpose of upper triangular R in same format as H
void qr(int N, double* H_array, double* R_array);

// Find inverse of hermitian NxN matrix H
// Same format for H as for qr function:
// H_00_re, H_10_re, H_10_im, H_11_re, H_20_re, H_20_im, H_21_re, H_21_im, ...
void inverse(int N, double* H_array, double* inverseH_array);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */