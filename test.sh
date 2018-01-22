#specify compilers if needed:
#export CXX=g++
#export CC=gcc

#compile qr library:
$CXX qr.cpp -shared -Wall -fPIC -Wno-int-in-bool-context -O3 -o libqr.so

#compile c example program
$CC -std=c99 -Wall -O3 example.c -L. -lqr -o example

#run c example program
LD_LIBRARY_PATH=. ./example
