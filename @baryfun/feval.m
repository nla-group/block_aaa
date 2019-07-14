function R = feval(obj, z)
%FEVAL    Evaluate BARYFUN at a scalar argument.
%
% Calling syntax: 
%   - R = feval(obj, z)  -- evaluate at scalar z

N = zeros(size(obj.Ck{1})); % numerator
D = zeros(size(obj.Dk{1})); % denominator

[val,ind] = min(abs(z-obj.zk)); % evaluation at support point
if val < 10*eps
    R = obj.Dk{ind}\obj.Ck{ind};
    return
end

for j = 1:length(obj.zk)
    N = N + obj.Ck{j}/(z-obj.zk(j));
    D = D + obj.Dk{j}/(z-obj.zk(j));
end
R = D\N;

end

