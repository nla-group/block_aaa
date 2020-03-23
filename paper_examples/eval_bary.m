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
    R = Dk{ind}\Ck{ind};
    return
end

for j = 1:length(zk)
    N = N + Ck{j}/(z-zk(j));
    D = D + Dk{j}/(z-zk(j));
end
R = D\N;

end

