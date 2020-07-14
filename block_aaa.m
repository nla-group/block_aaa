function [R,rmse,opts] = block_aaa(F,pts,opts)
%BLOCK_AAA   Block-AAA algorithm for discrete rational LS approximation.
%
% Block generalization of the AAA algorithm presented in 
% [Nakatsukasa/Sete/Trefethen, SIAM J. Sci. Comput. 40 (3), 2018]
% to m-by-n functions. Block-AAA is based on a generalized barycentric
% formula with matrix-valued weights, that is,
%
%    R(z) = (sum_i W_i/(z-z_i))\(sum_i W_i F(z_i)/(z-z_i))
%
% where the W_i are m-by-m matrices and z_i are scalar support points. 
% Note that the F(z_i) are m-by-n matrices and hence the numerator 
% effectively uses m-by-n matrices for interpolation (if W_i is nonsing.)
%
% Block-AAA is called as [R,rmse,out] = block_aaa(F,pts,opts)
%
% Inputs:  F    -- function handle @(z) to the m-by-n function, 
%                  or a cell array of length(pts) m-by-n matrices 
%          pts  -- vector of distinct sampling points in the complex plane
%          opts -- structure of additional options (default values):
%                  opts.tol    : target root mean squared error (1e-12)
%                  opts.maxit  : maximal number of iterations 
%                  opts.return : which baryfun to return, (best), last, all.
%
% Returns: R    -- rational interpolant as a BARYFUN object
%          rmse -- root mean squared error for each block-AAA iteration
%          out  -- additional outputs.
%
% Website: see https://github.com/nla-group/block_aaa for more details  
%

% handle input options
if nargin < 3, opts = struct(); end
if ~isfield(opts,'tol'), opts.tol = 1e-12; end
if ~isfield(opts,'maxit'), opts.maxit = 20; end
if ~isfield(opts,'svd'), opts.svd = 1; end % use SVD, otherwise EIG (not recomm.)
if ~isfield(opts,'chol'), opts.chol = 0; end % use QR+Cholesky update, otherwise full SVD
if ~isfield(opts,'return'), opts.return = 'best'; end % can be 'best', 'last', 'all'

pts = pts(:);
npts = length(pts);

% if F is function handle, evaluate it at all the pts
if iscell(F)
    FF = F;
else
    FF = cell(npts,1);
    for i = 1:npts
        FF{i} = F(pts(i));
    end
end
[m,n] = size(FF{1});

% compute norm ||F(z)|| at all pts
err = 0*pts;
for i = 1:npts
    err(i) = norm(FF{i},'fro').^2;
end
% find max and assign first support point there
[~,ind] = max(err);

lamind = 1:npts;              % indices of approximation points
lamind(lamind==ind) = [];
zk_ind = ind;
lam = pts(lamind);
zk = pts(zk_ind);

% M is the block Loewner matrix
% When updating is used, Qu and H are its QR factors such that M = Qu*H
M = zeros(0,(npts-1)*n);
Qu = zeros(0,(npts-1)*n); H = [];

if strcmp(opts.return,'all'), R = cell(1); else R = 0; end

for it = 1:opts.maxit+1
    
    ell = length(zk)-1; % degree
    p = npts-ell-1;
    
    if ~opts.chol   % compute SVD of full Loewner matrix at every iteration
        for i = 1:p % add new row to M
            M(ell*m+1:(ell+1)*m, (i-1)*n+1:i*n) = (FF{lamind(i)} - FF{zk_ind(ell+1)})/(lam(i) - zk(ell+1));
        end
        %[U,~,~] = svd(M,0);    % this is slow
        [~,~,U] = svd(M',0);    % tall skinny is faster!
        EF = U(:,end+1-m:end)';
    else % use QR+Chol update
        % add one more row to Qu (and M)
        for i = 1:p
            %M(ell*m+1:(ell+1)*m, (i-1)*n+1:i*n) = (FF{lamind(i)} - FF{zk_ind(ell+1)})/(lam(i) - zk(ell+1));
            Qu(ell*m+1:(ell+1)*m, (i-1)*n+1:i*n) = (FF{lamind(i)} - FF{zk_ind(ell+1)})/(lam(i) - zk(ell+1));
        end
        % orthogonalize m final rows of Qu (M = H*Qu)
        for r = ell*m+1:(ell+1)*m
            for i = 1:r-1
                H(r,i) = Qu(r,:)*Qu(i,:)';
                Qu(r,:) = Qu(r,:) - H(r,i)*Qu(i,:);
            end
            H(r,r) = norm(Qu(r,:));
            Qu(r,:) = Qu(r,:)/H(r,r);
        end
        [U,~,~] = svd(H,0);
        EF = U(:,end+1-m:end)';
    end
    
    % construct function handle to barycentric evaluation
    Ck = cell(ell+1,1); Dk = cell(ell+1,1);
    for j = 1:ell+1
        Dk{j} = EF(:,(j-1)*m+1:j*m);
        Ck{j} = Dk{j}*FF{zk_ind(j)};
    end
    % we don't use baryfun class here for performance
    Rbary = @(z) eval_bary(z,zk,Ck,Dk); 
    
    % evaluate RMSE
    err = 0*pts;
    for i = 1:npts
        err(i) = norm(Rbary(pts(i)) - FF{i},'fro').^2;
    end
    
    indnan = find(not(isfinite(err)));
    err(indnan) = 0;
    
    rmse(it) = sqrt(sum(err)/npts);
    
    if strcmp(opts.return,'best') && rmse(it) <= min(rmse)
        R = baryfun(zk,Ck,Dk);
    end
    if strcmp(opts.return,'last')
        R = baryfun(zk,Ck,Dk);
    end
    if strcmp(opts.return,'all')
        R{it} = baryfun(zk,Ck,Dk);
    end
    
    if rmse(it) < opts.tol || it == opts.maxit
        break
    end
    
    % not done: find max error and assign next support point
    [~,ind] = max(err);
    ind2 = find(lamind==ind);
    lamind(ind2) = [];
    zk_ind = [zk_ind , ind];
    lam = pts(lamind);
    zk = pts(zk_ind);
    
    if ~opts.chol
        M(:,1+(ind2-1)*n:ind2*n) = [];
    else
        Qc = Qu(:,1+(ind2-1)*n:ind2*n); % columns removed from Qu
        Qu(:,1+(ind2-1)*n:ind2*n) = [];
        % we still have:
        %disp(['before: ' num2str(norm(H*Qu - M)/norm(M))])
        % but Qu no longer orthonormal
        % update Qu
        Mat = (1+5*eps)*eye(it*m) - Qc*Qc';
        [S,flag] = chol(Mat);
        if flag
            disp(['CHOL warning, resorting to QR at iteration ' num2str(it)])
            [Q,R] = qr(Qu',0);
            H = H*R'; Qu = Q';
        else
            Qu = S'\Qu;
            H = H*S';
        end
        %disp(['after: ' num2str(norm(H*Qu - M)/norm(M))])
        %disp(['ortho: ' num2str(norm(Qu*Qu' - eye(it*m)))])
    end
    
end % for it
opts.zk = zk;
end % block_aaa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = eval_bary(z,zk,Ck,Dk)
%EVAL_BARY   Evaluate matrix-valued barycentric formula
%   
%      R(z) = inv(sum_k Dk/(z-zk))*(sum_k Ck/(z-zk))
%
% where z is the scalar evaluation point, zk is a vector of the distinct
% support points, and Ck and Dk are cell arrays with the corresponding
% numerator and denominator coefficients, respectively. 

N = zeros(size(Ck{1})); % numerator
D = zeros(size(Dk{1})); % denominator

[val,ind] = min(abs(z-zk)); % evaluation at support point
if val < 1e1*eps
    R = (Dk{ind})\Ck{ind};
    return
end

for j = 1:length(zk)
    N = N + Ck{j}/(z-zk(j));
    D = D + Dk{j}/(z-zk(j));
end
R = D\N;

end

