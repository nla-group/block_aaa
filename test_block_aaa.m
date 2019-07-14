clear, close all; clc

%% Example 1: approximate square 3x3 matrix F(z) on imaginary axis
load test_block_aaa
% provides FF, a cell array of 3x3 matrices (function values)
% and also pts, a vector of the sampling points 
% the data corresponds to ex2_Y.m in the Matrix Fitting Toolbox
% [B. Gustavsen and A. Semlyen, Rational approximation of frequency 
%  domain responses by Vector Fitting, IEEE Trans. Power Delivery 14(3),
%  pp. 1052--1061, 1999]

opts.tol = 1e-8;
opts.maxit = 21;
opts.return = 'best';
[R1,rmse1,out1] = block_aaa(FF,pts,opts); 
disp(['EXAMPLE 1 - best RMSE achieved: ' num2str(rmse(pts,R1,FF)) ])
figure, semilogy(0:length(rmse1)-1,rmse1), hold on
semilogy([0,opts.maxit-1],opts.tol*[1,1],'k--')
xlabel('degree'), ylabel('RMSE')
title('block-AAA convergence (Example 1, square F)')

%% Example 2: approximate rectangular 2x3 matrix F(z)
FF = @(z) [ 1/(z+1) 1/(z^2-5) z;
    1/(z^2+5+1i) (2+z^2)/(z^3 + 3*z^2 + 1) 7];
[R2,rmse2,out] = block_aaa(FF,pts,opts); 
disp(['EXAMPLE 2 - best RMSE achieved: ' num2str(rmse(pts,R2,FF)) ])
figure, semilogy(0:length(rmse2)-1,rmse2), hold on
semilogy([0,opts.maxit-1],opts.tol*[1,1],'k--')
xlabel('degree'), ylabel('RMSE')
title('block-AAA convergence (Example 2, rectangular F)')
