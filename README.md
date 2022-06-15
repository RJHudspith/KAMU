# KAMU

Here we have the numerical integration of the QED kernel to test against known results for the lepton loop.

The code includes hcubature.c and some related files from the cubature library \href{https://github.com/stevengj/cubature} this is licensed under GPL V2, which is in accord with our own license GPL V3.

The code requires a separate installation of the KQED library and then this code should be linked to it with the --with-KQED=<> configure option.

======= Some Examples =========

To reproduce the plots Fig.18 and Fig.19 from https://arxiv.org/pdf/2006.16224.pdf run:

   	./KINT 1 6E-3 MKERN 0.05 false

(as Method 1 is more peaked at small y this value is important) and

	./KINT 1 6E-3 MKERN 0.1 true

Don't set ymin or the tolerance much smaller than this or the code will have a bad time.