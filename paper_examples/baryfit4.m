function [Rbary,rmse,zk,Ck,Dk] = baryfit4(F,lam,iter)
%BARYFIT4    Surrogate-based AAA algorithm as described in
% S. Elsworth and S. Guettel, Conversions between barycentric, RKFUN, 
% and Newton representations of rational interpolants}, 
% Linear Algebra Appl., 576:246--257, 2019.
% 
% Applies the AAA algorithm to f(z)=a'*F(z)*b with random a,b.

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
a = rand(1,m); b = rand(n,1);
f = zeros(1,p);
for i = 1:p
    f(i) = a*FF{i}*b;
end

[r,pol,res,zer,zk,fk,wk,errvec] = aaa(f, lam, 0, iter+1);

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

