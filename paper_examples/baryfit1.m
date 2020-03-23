function [Rbary,rmse,zk,Ck,Dk] = baryfit1(F,lam,ell,opts)
%BARYFIT1    Non-interpolatory block fitting. Uses ell+1 random convex 
% combinations of pairs of sampling points lam as support points zk.
%
% block-NONFIT

if nargin < 4
    opts.svd = 1;
    opts.scale = 0;
    opts.plot = 0;
    opts.iter = 1;
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
[m,n] = size(FF{1});

% select ell+1 random support points
ind = randperm(p,2*(ell+1));
ind1 = ind(1:ell+1); ind2 = ind(ell+2:end);
alph = rand(ell+1,1); 
zk = alph.*lam(ind1) + (1-alph).*lam(ind2) + 5;

%zk = lam(randperm(p,ell+1));



% form Cauchy matrix
[Lam,Zk] = meshgrid(lam,zk);
C = 1./(Lam - Zk); 

Color = { 'g','b','r','c','m' };
for it = 1:opts.iter
    M1 = kron(C,-eye(n)); 
    M2 = zeros((ell+1)*m,p*n);
    for i = 1:p
        M2(:,(i-1)*n+1:i*n) = kron(C(:,i),FF{i});
    end
    M = [ M1 ; M2 ];
    
    if opts.scale % column scaling of M
        d = diag(1./sqrt(sum(abs(M).^2,1)));
        M = M*d;
    end
    
    % compute left singular block
    if opts.svd
        [U,~,~] = svd(M,0);
        EF = U(:,end+1-m:end)';
    else
        [U,D] = eig(M*M'); 
        [~,ind] = sort(diag(D));
        EF = U(:,ind(1:m))';
    end
    
    % function handle to barycentric evaluation
    Ck = cell(ell+1,1); Dk = cell(ell+1,1);
    for j = 1:ell+1
        Ck{j} = EF(:,(j-1)*n+1:j*n); 
        Dk{j} = EF(:,(ell+1)*n+(j-1)*m+1:(ell+1)*n+j*m);
    end 
    Rbary = @(z) eval_bary(z,zk,Ck,Dk);

    err = 0*lam;
    for i = 1:p
        err(i) = norm(Rbary(lam(i)) - FF{i},'fro').^2;
    end
    if opts.plot
        figure
        semilogy(1:p,sqrt(err),'-','Color',Color{1+mod(it,5)})
        xlabel('# sampling point')
        ylabel('absolute error in Frobenius norm')
    end

end % iter

rmse = sqrt(sum(err)/p);
return


%% form rkfunb (ignore)
H = full(kron(spdiags([[zk(2:end);nan],-zk],-1:0,ell+1,ell),speye(s)));
K = full(kron(spdiags([ones(ell+1,1),-ones(ell+1,1)],-1:0,ell+1,ell),speye(s)));
C = mat2cell(EF(:,1:(ell+1)*s),s,repmat(s,ell+1,1)).';
D = mat2cell(EF(:,1+(ell+1)*s:end),s,repmat(s,ell+1,1)).';
C = cell2mat(C); D = cell2mat(D);
[QQ,RR] = qr(D); QQ(:,1:s) = QQ(:,1:s)*RR(1:s,1:s);
H = QQ\H; K = QQ\K; C = QQ\C;
[~,~,QQ,ZZ] = qz(H(s+1:end,:),K(s+1:end,:));
QQ = blkdiag(eye(s),QQ);
K = triu(QQ*K*ZZ,-s); H = triu(QQ*H*ZZ,-s); C = QQ*C;
C = mat2cell(C,repmat(s,1,ell+1),s);
R = rkfunb(K, H, C, 0);

if opts.plot
    figure
    plot(real(lam),imag(lam),'k.'), hold on
    plot(real(zk),imag(zk),'ko')
    plot(roots(R),'bx'), plot(poles(R),'r+')
    legend('sampling points','support points','roots','poles')
end
