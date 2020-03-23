Instructions file for using the Matlab codes associated to the report

"Algorithms for the rational approximation of matrix-valued functions"

by Ion Victor Gosea and Stefan Guettel

available online at 

https://arxiv.org/pdf/2003.06410.pdf

-------------------------------------------------------------------------------------------------------------------

Dependent pacakges

1) Chebfun package

The Chebfun package (containing MATLAB implementations of the AAA 
algorithm and many more) is available at

https://www.chebfun.org/

Mored details at:

T. A. Driscoll, N. Hale, and L. N. Trefethen, Chebfun Guide. 
https://www.chebfun.org/docs/guide/, 2014.


2) Rational Krylov Toolbox

A MATLAB implementation of RKFIT is provided at 

http://www.rktoolbox.org/

More details at:

M. Berljafa, S. Elsworth, and S. Guettel, A Rational Krylov Toolbox 
for MATLAB, MIMS EPrint 2014.56, Manchester Institute for Mathematical
Sciences, The University of Manchester, UK, 2014.

3) Matrix Fitting: 

The extension of vector fitting to matrix-valued data was implemented 
in the MATLAB Matrix Fitting Toolbox available at

https://www.sintef.no/projectweb/vectorfitting/downloads/matrix-fitting-toolbox/

One also requires the mat file "ex2_Y.mat" from the toolbox above 
(to be used in example2.m)

More details at:

B. Gustavsen, Improving the pole relocating properties of vector fitting, 
IEEE Trans. Power Delivery, 21 (2006), pp. 1587-1592.

4) Block-AAA:

A MATLAB implementation of the block-AAA method described in the arxiv
report can be downloaded from

https://github.com/nla-group/block_aaa

5) Benchmark Examples for Model Reduction

The linear models needed for running the demo files example3.m, 
example4.m and example6.m (corresponding to the experiments in 
Section 5.3, 5.4 and 5.6 in the arxiv report). 
The mat files "CDplayer.mat" and "iss.mat" are available at 

http://slicot.org/20-site/126-benchmark-examples-for-model-reduction

----------------------------------------------------------------------------------------------------------------------------------------------------

Main DEMO Files

1) example1.m 

- the test cases are symmetric/nonsymmetric toy examples (Section 5.1 in 
the arxiv report).

2) example2.m

- the test case is a collection of 300 samples of 3-by-3 admittance 
matrices for a three-conductor overhead line from the Matrix Fitting 
Toolbox (Section 5.2 in the arxiv report).

3) example3.m

- the test case is the CD player model from the SLICOT benchmark 
collection  (Section 5.3 in the arxiv report).

4) example4.m

- the test case is the ISS model from the SLICOT benchmark collection  
(Section 5.4 in the arxiv report).

5) example5.m

- the test case is the buckling plate model  (Section 5.5 in the arxiv 
report) from

N. J. Higham, G. M. Negri Porzio, and F. Tisseur, An updated set of 
nonlinear eigenvalue problems, Tech. Rep. MIMS EPrint: 2019.5, Manchester,
United Kingdom, 2019. http://eprints.maths.manchester.ac.uk/.

6) example6.m

- the test case is a scalar example with noise (Section 5.6 in the 
arxiv report).

7) example7.m

- the test case is the ISS model also used in example4.m, but this 
time noise is added to the measurements (Section 5.7 in the arxiv report).


