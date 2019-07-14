function D = nonlinear_eig(C,z)

[m,n] = size(C{1});
if m~=n
    error('nonlinear eigenvalue/root finding works only for square BARYFUNs');
end

ell = length(z)-1;

if ell <= 0
    D = [];
    return
end

% find weights w_j of barycentric polynomial interpolant
w = zeros(ell+1,1);
for i=1:ell+1
    x = z - z(i);
    x(i) = [];
    w(i) = 1/prod(x);
    C{i} = C{i}/w(i); % rescale C's
end
%w = ones(1,ell+1);

theta = w(1:end-1)./w(2:end);

% Split L0 into top and bottom matrices
A = zeros(ell-1,ell);
for i=1:ell-1
    A(i,i) = z(i);
    A(i,i+1) = -z(i+2)*theta(i);
end
L0bott = kron(A,eye(m,m)); % Form bottom matrix in L0

% Top matrix in L0
L0top = [];
for i=1:ell-1
   L0top = [L0top, z(i+1)*C{i}]; 
end
L0top = [L0top, z(ell+1)*C{ell} + z(ell)*(1/theta(ell))*C{ell+1}];
% Concatenate to create L0
L0 = [L0top; L0bott];
 
% Similar computation for L1
A = zeros(ell-1,ell);
for i=1:ell-1
    A(i,i) = 1;
    A(i,i+1) = -theta(i);
end
L1bott = kron(A,eye(m,m));

L1top = [];
for i=1:ell-1
   L1top = [L1top, C{i}];
end

L1top = [L1top, C{ell} + (1/theta(ell))*C{ell+1}];
% Concatenate to create L1
L1 = [L1top; L1bott];

% Find generalized eigvalues of (L0,L1)
D = sort(eig(L0,L1));


