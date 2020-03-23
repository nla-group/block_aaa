function [R,rmse] = baryfit6(F,lam,xi,ops)
% RKFIT wrapper to provide F as a matrix with varying functions for the entries.
% A matrix of rational functions 'ratfun' is returned, all sharing the same
% common denominator, with poles xi.
%
% if xi = [], then ops.iter incremental RKFIT iteration will be run,
% adding one pole per iteration.


% Initial values
npts = length(lam); % length of initial support points
A = spdiags(lam(:),0,npts,npts); % Populate matrix A
b = ones(npts,1); % Initial weighting of ones for b

if nargin<4
    ops.iter = 3;
end

if isempty(xi)    % fixed by SG on 14/02/2020 
    ops.iter = 0;
end

if iscell(F) % Check if F is already a cell array of populated data
    FF = F;
else  % If not, populate the cell with data lam using matrix F
    FF = cell(1,npts);
    for i=1:npts
        FF{i} = F(lam(i));
    end
end

    [m,n] = size(FF{1});
    Flam = zeros(npts,m*n);
    for i=1:npts
        Flam(i,:) = reshape(FF{i},1,m*n);
    end
    FFF = cell(1,m*n);
    for i=1:m*n
        FFF{1,i} = spdiags(Flam(:,i),0,npts,npts);
    end

param.maxit = ops.iter;
param.reduction = 0;
param.tol = 0;
param.return = 'best';
param.reorth = 0;
[xi_out,ratfun,misfit,out] = rkfit(FFF,A,b,xi,param);

% 
if length(ratfun) == 1 % If F is a single input function just use ratfun
    R = ratfun;
else % Evaluation of the function ratfun in matrix form
    %R = @(z) reshape(cellfun(@(fun) feval_qr(fun,z),ratfun),m,n);
    %% more efficient eval of multiple ratfuns
    K = ratfun{1}.K;
    H = ratfun{1}.H;
    coeffs = zeros(length(ratfun{1}.coeffs),length(ratfun));
    for j = 1:length(ratfun)
        coeffs(:,j) = ratfun{j}.coeffs;
    end
    R = @(z) reshape(eval_ratfuncell(z,K,H,coeffs),m,n);
end



%%
err = zeros(1,npts);
for i=1:npts
    err(i) = norm(R(lam(i))-FF{i},'fro')^2;
end
% Calculate the RMSE
rmse = sqrt(sum(err)/npts);


function v = eval_ratfuncell(z,K,H,coeffs)
% used to avoid multiple QR calls when evaluating cell array of ratfuns
n = null((z*K - H).').'; 
n = n/n(1);
v = n*coeffs;


