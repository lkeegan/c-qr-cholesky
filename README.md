# c-qr-cholesky [![Build Status](https://travis-ci.org/lkeegan/c-qr-cholesky.svg?branch=master)](https://travis-ci.org/lkeegan/c-qr-cholesky)
C interface to QR and Cholesky matrix decomposition routines using the [Eigen](http://eigen.tuxfamily.org) C++ template library for matrix operations (included)

Based on: https://stackoverflow.com/questions/26889142/using-eigen-in-a-c-project

To compile the library:

```
g++ qr.cpp -shared -Wall -fPIC -O3 -o libqr.so
```

To compile the example c program that calls the library:

```
gcc -std=c99 -Wall -O3 -lm example.c -L. -lqr -o example
```

To run the example c program:

```
LD_LIBRARY_PATH=. ./example
```

Alternatively all three steps are included in the file ./test.sh