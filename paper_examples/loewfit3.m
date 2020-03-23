function [R,rmse,zk,pol,zer] = loewfit3(F,lam,ell,opts)
%
%
% This function implements the barycentric Loewner rational approx algorithm
% for the MIMO (multi-input multi-output) case
%
% It constructs a rational function R (in function handle format) that
% approximates the data pair (F(lam),lam)
% 
%% INPUT arguments
%
% F = function to be approximated (or a vector of samples values)
% lam = sampling points (on which the function F is approximated)
% ell = order of the rational fct. that approximates F: (r-1,r) to be exact 
% opts = options (partition type, representation type or plotting (yes or no)) 
%
%
% opts.partition = 1 => interlacing points (alternate partition);
% opts.partition = 0 => separated points (half-half partition);
%
% opts.representation = 1 => Loewner rational function (R) in barycentric format;
% opts.representation = 0 => Loewner rational function (R) in transfer function format;
%
% opts.directions = 1 => tangential directions chosen randomly
% opts.directions = 0 => tangential directions chosen as ones vectors
%
%% OUTPUT arguments
%
% R    = rational interpolant computed using the Loewner framework
% rmse = the total RMSE of R at all the sampling points
% zk   = the vector of support points
% Ck   = cell array containing the numerator matrices
% Dk   = cell array containing the denominator matrices
% Pk   = cell array containing the support point matrices
% pol  = the poles of R
% zer  = the zeros of R

% Ion Victor Gosea, MPI Magdeburg
%
% Last revision 22 July 2019


if nargin < 4
    opts.partition = 1;
    opts.plot = 0;
    opts.representation = 0;
    opts.directions = 1;
    opts.svd = 1;
end


% try two different partions of the interpolation points (alternate(1) or
% half-half(0))
if(opts.partition)
            left_ind  = 1:2:length(lam);
            right_ind = 2:2:length(lam);
else      
            left_ind  = 1:length(lam)/2;
            right_ind = (length(lam)/2+1):length(lam);
end

lam = lam(:);
p = length(lam);

% if F is function handle, evaluate it at all the lam's
if iscell(F)
    FF = F;
else
    FF = cell(p,1);
    for i = 1:p
        FF{i} = F(lam(i));
    end
end

m = size(FF{1},1);
n = size(FF{1},2);

k = p/2;

V = zeros(k,n);
W = zeros(m,k);

if(opts.directions)
    LM = randn(k,n);
    RM = randn(n,k);
else
    LM = ones(k,n);
    RM = ones(n,k); 
end

for jj = 1:k
    V(jj,1:n) = LM(jj,:)*FF{left_ind(jj)};
    W(1:m,jj) = FF{right_ind(jj)}*RM(:,jj);
end


% repeat the interpolation points 
Leftp  = lam(left_ind).';
Rightp = lam(right_ind);

% diagonal matrices with the left and right points on the main diagonal
LeftM = diag(Leftp);
RightM = diag(Rightp);

% Put together the Loewner matrix
%LL  = (V*RightM-LeftM*W)./(V*RM-LM*W);

LL = sylvester(LeftM,-RightM,V*RM-LM*W);

% Put together the shifted Loewner matrix
LLs = sylvester(LeftM,-RightM,LeftM*V*RM-LM*W*RightM);

if(opts.svd)

    % Perform a SVD of the Loewner matrix
    [Xs,SV,Ys]=svd(LL);
         
    % collect the normalized singular values in vector format
    sv = diag(SV); sv = sv/sv(1);

    if(opts.plot)
       semilogy(sv);
       title('Singular value decay of the Loewner matrix');
    end

    % Put together the projection matrices (the new way)
    X = Xs(:,1:ell)  ;
    Y = Ys(:,1:ell)  ;
    %S = SV(1:ell,1:ell);

else
    
    [Xs1,SV1,Ys1]=svd([LL LLs],'econ');
    [Xs2,SV2,Ys2]=svd([LL;LLs],0);
    
    X = Xs1(:,1:ell)  ;
    Y = Ys2(:,1:ell)  ;
end

% Projected quantities
Er = - X'*LL*Y;
Ar = - X'*LLs*Y;
Br = X'*V;
Cr = W*Y;

% the rational interpolant
R = @(x) Cr*((x*Er-Ar)\Br);

% the vector of support points
zk = eig(Y'*RightM*Y);

% Compute the poles of the rational interpolant
pol = eig(Ar);
 
% Compute the zeros of the rational interpolant
zer = eig([Ar Br;Cr zeros(m,n)],blkdiag(eye(ell),zeros(m,n)));

for ii = 1:p
    err(ii) = norm(R(lam(ii)) - FF{ii},'fro').^2;
end

rmse = sqrt(sum(err)/p);


end