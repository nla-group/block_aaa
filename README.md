# Block-AAA
## A method for discrete rational least squares approximation

This repository provides a MATLAB implementation of the block-AAA method
described in [3]. This method is a block generalization of the AAA 
algorithm [7] to m-by-n matrix-valued functions F(z) of a scalar argument. 
Given a function F(z) at a number of sampling points in the complex plane, 
the block-AAA method aims to compute a low-degree rational interpolant R(z) 
such that R(z)~F(z), measured in the root mean squared error over the 
sampling points.

Block-AAA is based on a generalized barycentric formula with 
*matrix-valued weights,* that is,

   R(z) = inv(sum_i W_i/(z-z_i))*(sum_i W_i F(z_i)/(z-z_i))

where the W_i are m-by-m matrices and the z_i are scalar support points. 
Note that the F(z_i) are m-by-n matrices and hence the numerator of R(z) 
effectively uses m-by-n matrices. If W_i is nonsingular and the z_i
are all distinct (which we assume), then R(z_i) = F(z_i), i.e., 
R(z) interpolates F(z) at the support points z_i. 
In block-AAA, the support points z_i are chosen one after the other in 
a greedy fashion and the weights W_i are determined by linearized least 
squares approximation, as in the AAA algorithm [7]. 

Other matrix-valued extensions of AAA have been described in [2, 5, 6]. 
These methods produce rational interpolants with a *scalar denominator
polynomial*, whilst a block-AAA interpolant has a *matrix polynomial* in its 
denominator. Hence, a block-AAA interpolant of degree d can have up to 
m*d finite poles, whereas the scalar-denominator interpolants in [2, 5, 6] 
of the same degree have at most d finite poles. The additional flexibility 
of the block-AAA interpolant can yield reduced approximation degrees 
for the same accuracy. On the other hand, such interpolants are more 
difficult to fit and block-AAA can exhibit erratic convergence behaviour. 
Further, many spurious poles and nonlinear eigenvalues might appear. 

A comprehensive discussion and comparison of many methods for discrete 
rational least squares approximation, including the direct use of the 
Loewner approach [1] and vector fitting [4], will be subject of [3].


## Example

`block_aaa` requires two essential inputs, a function handle or cell array of function values `FF`, and a vector of sampling points `pts`. Further options can be provided in a structure `opts`.

```matlab
>> load test_block_aaa
>> opts.tol = 1e-8;
>> opts.maxit = 21;
>> [R1,rmse1,out1] = block_aaa(FF,pts,opts); 
>> rmse(pts, FF, R1)

ans =

   6.0822e-09
```

See also the file `test_block_aaa.m` in the main directory of the repository.


## License

This project is licensed under the MIT License.


## References

[1] A. C. Antoulas and B. D. Q. Anderson: On the scalar rational 
    interpolation problem, IMA Journal of Mathematical Control and 
    Information 3, pp. 61–-88, 1986.

[2] S. Elsworth and S. Güttel: Conversions between barycentric, RKFUN, 
    and Newton representations of rational interpolants, 
    Linear Algebra and its Applications 576, pp. 246--257, 2018.

[3] I. V. Gosea and S. Güttel: Algorithms for the rational approximation 
    of matrix-valued functions, accepted for publication in SIAM Journal 
    on Scientific Computing, 2021. [https://arxiv.org/abs/2003.06410]

[4] B. Gustavsen and A. Semlyen: Rational approximation of frequency 
    domain responses by vector fitting, IEEE Transactions on Power 
    Delivery 14(3), pp. 1052--1061, 1999.

[5] A. Hochman. FastAAA: A fast rational-function fitter. 
    In: 2017 IEEE 26th Conference on Electrical Performance of Electronic 
    Packaging and Systems (EPEPS), pp. 1--3, 2017.

[6] P. Lietaert, J. Pérez, B. Vandereycken, and K. Meerbergen: Automatic 
    rational approximation and linearization of nonlinear eigenvalue 
    problems, arXiv preprint 1801.08622, 2018.
    [https://arxiv.org/pdf/1801.08622.pdf]

[7] Y. Nakatsukasa, O. Sète, and L. N. Trefethen: The AAA algorithm for 
    rational approximation, SIAM Journal on Scientific Computing 40(3), 
    pp. A1494--A1522, 2018.
