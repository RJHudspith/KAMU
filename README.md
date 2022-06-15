# KAMU

Here we have the numerical integration of the QED kernel to test against known results for the lepton loop.

The code includes hcubature.c and some related files from the cubature library \href{https://github.com/stevengj/cubature} this is licensed under GPL V2, which is in accord with our own license GPL V3.

The code requires a separate installation of the KQED library and then this code should be linked to it with the --with-KQED=<> configure option.
