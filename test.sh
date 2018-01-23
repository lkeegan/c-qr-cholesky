#specify compilers if needed:
#export CXX=clang++
#export CC=clang

#compile qr library:
$CXX qr.cpp -shared -Wall -fPIC -Wno-int-in-bool-context -O3 -o libqr.so

#compile c example program
$CC -std=c99 -Wall -O3 -lm example.c -L. -lqr -o example

#run c example program
LD_LIBRARY_PATH=. ./example
