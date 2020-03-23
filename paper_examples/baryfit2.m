function [Rbary,rmse,zk,Ck,Dk] = baryfit2(F,pts,iter,opts)
%BARYFIT2   Block AAA algorithm (interpolatory).

if nargin < 4
    opts.scale = 0; % deprecated
    opts.plot = 0;
    opts.chol = 0; % use QR + Cholesky update, otherwise full SVD
    opts.svd = 1; % use SVD, otherwise eig (not recommended)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pts = pts(:);
npts = length(pts);

% if F is function handle, evaluate it at all the lam's
if iscell(F)
    FF = F;
else
    FF = cell(npts,1);
    for i = 1:npts
        FF{i} = F(pts(i));
    end
end
[m,n] = size(FF{1});

Rbary = @(z) zeros(m,n);     % initial approximant

% compute error at approximation points
err = 0*pts;
for i = 1:npts
    err(i) = norm(Rbary(pts(i)) - FF{i},'fro').^2;
end
% find max and assign first support point
[~,ind] = max(err);

lamind = 1:npts;              % indices of approximation points
lamind(lamind==ind) = [];
zk_ind = ind;
lam = pts(lamind);
zk = pts(zk_ind);

Color = { 'b','r','c','m','g' };

% for updating the Loewner matrix (and its QR factor)
M = zeros(0,(npts-1)*n);
Qu = zeros(0,(npts-1)*n); H = [];

for it = 1:iter+1
    
    ell = length(zk)-1; % degree
    p = npts-ell-1;
    
    if ~opts.chol  % compute SVD of Loewner matrix at every iteration
        for i = 1:p  % just add new row to M
            M(ell*m+1:(ell+1)*m, (i-1)*n+1:i*n) = (FF{lamind(i)} - FF{zk_ind(ell+1)})/(lam(i) - zk(ell+1));
        end
        if opts.svd % use SVD (instead of eig)
            %[U,~,~] = svd(M,0);
            [~,~,U] = svd(M',0);  % tall skinny is faster!
            EF = U(:,end+1-m:end)';
        else
            [U,D] = eig(M*M'); 
            [~,ind] = sort(diag(D));
            EF = U(:,ind(1:m))';
        end
        
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
    
    % function handle to barycentric evaluation
    Ck = cell(ell+1,1); Dk = cell(ell+1,1);
    for j = 1:ell+1
        Dk{j} = EF(:,(j-1)*m+1:j*m);
        Ck{j} = Dk{j}*FF{zk_ind(j)};
    end 
    
    Rbary = @(z) eval_bary(z,zk,Ck,Dk);
    
    err = 0*pts;
    for i = 1:npts
        err(i) = norm(Rbary(pts(i)) - FF{i},'fro').^2;
    end
    rmse(it) = sqrt(sum(err)/npts);
   
    
    if opts.plot
        figure(99)
        semilogy(1:npts,err,'-','Color',Color{1}), hold on
        semilogy(zk_ind,err(zk_ind),'o','Color',Color{1}), 
        xlabel('# sampling point'), ylabel('abs error in Frobenius norm')
        legend('error','support point'), shg
    end
    
    % if not done: find max error and assign next support point
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
        
        Mat = (1+0*eps)*eye(it*m) - Qc*Qc';
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
    
end % iter
zk = zk(1:iter+1);

if opts.plot
    figure
    semilogy(rmse)
    xlabel('Block AAA iteration'), ylabel('RMSE')
end