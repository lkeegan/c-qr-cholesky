#compile qr library:
g++ qr.cpp -shared -Wall -fPIC -Wno-int-in-bool-context -O3 -o libqr.so

#compile c example program
gcc example.c -L. -lqr -o example

#run c example program
LD_LIBRARY_PATH=. ./example
