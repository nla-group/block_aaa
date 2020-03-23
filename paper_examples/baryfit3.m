function [Rbary,rmse,zk,Ck,Dk] = baryfit3(F,lam,iter)
%BARYFIT3    Set-valued AAA algorithm by Lietart et al
% See https://arxiv.org/pdf/1801.08622.pdf 

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
Frow = zeros(p,m*n);
for i = 1:p
    Frow(i,:) = reshape(FF{i},1,m*n);
end

[Rrow, POL, RES, ZER, zk, fk, wk, err] = aaa_sv(Frow, lam, 'mmax', iter+1, 'tol', 0);

Ck = cell(length(wk),1); Dk = cell(length(wk),1);
for j = 1:length(wk)
    ind = find(lam == zk(j),1,'first');
    Ck{j} = wk(j)*FF{ind};
    Dk{j} = wk(j);
end
Rbary = @(z) eval_bary(z,zk,Ck,Dk);


% compute rmse
for i = 1:p
   err(i) = norm(FF{i}-Rbary(lam(i)),'fro').^2; 
end
rmse = sqrt(sum(err)/p);


end

