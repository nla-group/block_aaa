function [r, pol, res, zer, z, f, w, errvec] = aaa_sv(F, varargin)
%AAA   Computes a AAA rational approximation where F can have multiple
%outputs
%   R = AAA(F, Z) computes the AAA rational approximant R (function handle) to
%   data F on the set of sample points Z.  F may be given by its values at Z,
%   or as a function handle or a chebfun.
%
%   [R, POL, RES, ZER] = AAA(F, Z) returns vectors of poles POL,
%   residues RES, and zeros ZER of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ] = AAA(F, Z) also returns the vectors
%   of support points ZJ, function values FJ, and weights WJ of the
%   barycentric representation of R.
%
%   [R, POL, RES, ZER, ZJ, FJ, WJ, ERRVEC] = AAA(F, Z) also returns the
%   vector of errors ||f-r||_infty in successive iteration steps of AAA.
%
%   R = AAA(F, Z, NAME, VALUE) sets the following parameters:
%   - 'tol', TOL: relative tolerance (default TOL = 1e-13),
%   - 'mmax', MMAX: maximal number of terms in the barycentric representation
%       (default MMAX = 100).
%   - 'dom', DOM: domain (default DOM = [-1, 1]). No effect if Z is provided.
%   - 'cleanup', 'off' or 0: turns off automatic removal of numerical Froissart
%       doublets
%
%   One can also execute R = AAA(F), with no specification of a set Z.
%   This is equivalent to defining Z = LINSPACE(DOM(1), DOM(2), LENGTH(F)) if F
%   is a vector (by default DOM = [-1, 1]).
%   If F is a function handle or a chebfun, AAA attempts to resolve F on its
%   domain DOM.  By default, DOM = [-1, 1] for a function handle, and
%   DOM = F.DOMAIN([1, END]) for a chebfun.
%

% Parse inputs:
[F, Z, M, nF, dom, tol, mmax, cleanup_flag, needZ, mmax_flag] = parseInputs(F, varargin{:});

if ( needZ )
    % Z was not provided.  Try to resolve F on its domain.
    [r, pol, res, zer, zj, fj, wj, errvec] = ...
        aaa_autoZ(F, dom, tol, mmax, cleanup_flag, mmax_flag);
    return
end

% Scale the functions
normF = max(abs(F),[],1);
F = bsxfun(@rdivide,F,normF);

% Left scaling matrix:
SF = spdiags(F, 0:-M:-M*(nF-1), M*nF, M);

% Initialize values
F = F(:);
R = mean(F);
errvec = zeros(mmax,1);
z = zeros(mmax,1);
f = zeros(mmax,nF);
ind = zeros(mmax,nF);
H = zeros(mmax,mmax-1);
S = zeros(mmax,mmax-1);

Q = zeros(M*nF,0);
C = zeros(M,0);

% AAA iteration:
for m = 1:mmax
    [errvec(m),loc] = max(abs(F-R));               % Select next support point where error is largest
    if ( errvec(m) <= tol )
        m = m-1;
        break
    end
    
    loc = mod(loc,M);
    ind(m,:) =  loc + (M*(loc==0):M:(nF-1+(loc==0))*M);  % Get indices of the z_i
    z(m) = Z(ind(m,1));                           % Add interpolation point
    f(m,:) = F(ind(m,:));                         % Add function values
    C(:,end+1) = 1./(Z - z(m));                   % Get column of the Cauchy matrix.
    C(ind(1:m,1),m) = 0;                          % Set the selected elements to 0
    
    v = C(:,m)*f(m,:);                            % Compute the next vector of the basis.
    v = SF*C(:,m)-v(:);
    
    % Update H and S to compensate for the removal of the rows
    q = Q(ind(m,:),1:m-1);
    q = q*S(1:m-1,1:m-1);
    Si = chol(eye(m-1,m-1)-q'*q);
    H(1:m-1,1:m-1) = Si*H(1:m-1,1:m-1);
    S(1:m-1,1:m-1) = S(1:m-1,1:m-1)/Si;
    S(m,m) = 1;
    Q(ind(1:m,:),:) = 0;
    
    nv = norm(v);
    H(1:m-1,m) = Q'*v;
    H(1:m-1,m) = S(1:m-1,1:m-1)'*H(1:m-1,m);
    HH = S(1:m-1,1:m-1)*H(1:m-1,m);
    v = v-(Q*HH);
    H(m,m) = norm(v);
    % Reorthoganlization is necessary for higher precision
    it = 0;
    while (it < 3) && (H(m,m) < 1/sqrt(2)*nv)
        h_new = S(1:m-1,1:m-1)'*(Q'*v);
        v = v - Q*(S(1:m-1,1:m-1)*h_new);
        H(1:m-1,m) = H(1:m-1,m) + h_new;
        nv = H(m,m);
        H(m,m) = norm(v);
        it = it+1;
    end
    v = v/H(m,m);
    
    % Add v
    Q(:,end+1) = v;
    
    % Solve small least squares problem with H
    [~,~,V] = svd(H(1:m,1:m));
    w = V(:,end);
    
    % Get the rational approximation
    N = C*bsxfun(@times,w,f(1:m,:));       % Numerator
    D = C*bsxfun(@times,w,ones(m,nF));     % Denominator
    R = N(:)./D(:);
    R(ind(1:m,:)) = F(ind(1:m,:));
end

f = f(1:m,:);
w = w(1:m);
z = z(1:m);

% Scale function values back
f = bsxfun(@times,f,normF);

% Note: When M == 2, one weight is zero and r is constant.
% To obtain a good approximation, interpolate in both sample points.
if ( M == 2 )
    z = Z;
    f = F;
    w = [1; -1];       % Only pole at infinity.
    w = w/norm(w);   % Impose norm(w) = 1 for consistency.
    errvec(2) = 0;
end

% Remove support points with zero weight:
I = find(w == 0);
z(I) = [];
w(I) = [];
f(I,:) = [];

% Construct function handle:
r = @(zz) reval(zz, z, f, w);

% Compute poles, residues and zeros:
[pol, res, zer] = prz(r, z, f, w);

% if ( cleanup_flag )
%     % Remove Froissart doublets:
%     [r, pol, res, zer, zj, fj, wj] = cleanup(r, pol, res, zer, zj, fj, wj, Z, F);
% end

end % of AAA()



%% parse Inputs:

function [F, Z, M, nF, dom, tol, mmax, cleanup_flag, needZ, mmax_flag] = ...
    parseInputs(F, varargin)
% Input parsing for AAA.

% Check if F is empty:
if ( isempty(F) )
    error('CHEBFUN:aaa:emptyF', 'No function given.')
end

% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    if ( isempty(Z) )
        error('CHEBFUN:aaa:emptyZ', ...
            'If sample set is provided, it must be nonempty.')
    end
    varargin(1) = [];
end

% Set defaults for other parameters:
tol = 1e-13;        % Relative tolerance.
mmax = 100;         % Maximum number of terms.
% Domain:
if ( isa(F, 'chebfun') )
    dom = F.domain([1, end]);
else
    dom = [-1, 1];
end
cleanup_flag = 1;   % Cleanup on.
mmax_flag = 0;

% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            mmax = varargin{2};
            mmax_flag = 1;
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'dom', 3) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 2]) )
            dom = varargin{2};
        end
        varargin([1, 2]) = [];
        if ( isa(F, 'chebfun') )
            if ( ~isequal(dom, F.domain([1, end])) )
                warning('CHEBFUN:aaa:dom', ...
                    ['Given domain does not match the domain of the chebfun.\n', ...
                    'Results may be inaccurate.'])
            end
        end
        
    elseif ( strncmpi(varargin{1}, 'cleanup', 7) )
        if ( strncmpi(varargin{2}, 'off', 3) || ( varargin{2} == 0 ) )
            cleanup_flag = 0;
        end
        varargin([1, 2]) = [];
        
    else
        error('CHEBFUN:aaa:UnknownArg', 'Argument unknown.')
    end
end


% Deal with Z and F:
if ( ~exist('Z', 'var') && isfloat(F) )
    % F is given as data values, pick same number of sample points:
    Z = linspace(dom(1), dom(2), size(F,1)).';
end

if ( exist('Z', 'var') )
    % Z is given:
    needZ = 0;
    
    % Work with column vector:
    Z = Z(:);
    M = length(Z);
    
    % Function values:
    if ( isa(F, 'function_handle') || isa(F, 'chebfun') )
        % Sample F on Z:
        F = F(Z);
    elseif ( isnumeric(F) )
        % Work with column vector and check that it has correct length.
        if ( size(F,1) ~= M )
            error('CHEBFUN:aaa:lengthFZ', ...
                'Inputs F and Z must have the same length.')
        end
    else
        error('CHEBFUN:aaa:UnknownF', 'Input for F not recognized.')
    end
    
else
    % Z was not given.  Set flag that Z needs to be determined.
    % Also set Z and M since they are needed as output.
    needZ = 1;
    Z = [];
    M = length(Z);
end
nF = size(F,2);

end % End of PARSEINPUT().


%% Evaluate rational function in barycentric form.

function r = reval(zz, zj, fj, wj)
% Evaluate rational function in barycentric form.
l = length(zz);
zv = zz(:);                             % vectorize zz if necessary
CC = 1./bsxfun(@minus, zv, zj.');       % Cauchy matrix
r = (CC*bsxfun(@times,wj,fj))./(CC*wj);             % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv),:) = kron(ones(sum(isinf(zv)),1),sum(bsxfun(@times,wj,fj),1)./sum(wj));

% Deal with NaN:
ii = find(isnan(r));
ii = [mod(ii(:),length(zz)),floor(ii(:)/l)+1];
for jj = 1:size(ii,1)
    if ( isnan(zv(ii(jj,1))) || ~any(zv(ii(jj,1)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj,1),ii(jj,2)) = fj(zv(ii(jj,1)) == zj,ii(jj,2));
    end
end

% Reshape to input format:
% r = reshape(r, length(zz),size(fj,2));

end % End of REVAL().


%% Compute poles, residues and zeros.

function [pol, res, zer] = prz(r, zj, fj, wj)
% Compute poles, residues, and zeros of rational function in barycentric form.
m = length(wj);

% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

% Compute residues via discretized Cauchy integral:
% dz = 1e-5*exp(2i*pi*(1:4)/4);
% res = r(bsxfun(@plus, pol, dz))*dz.'/4;
res = 0;

% Compute zeros via generalized eigenvalue problem:
for it = 1:size(fj,2)
    E = [0 (wj.*fj(:,it)).'; ones(m, 1) diag(zj)];
    zer{it} = eig(E, B);
    % Remove zeros of numerator at infinity:
    zer{it} = zer{it}(~isinf(zer{it}));
end
end % End of PRZ().


%% Cleanup

function [r, pol, res, zer, z, f, w] = cleanup(r, pol, res, zer, z, f, w, Z, F)
% Remove spurious pole-zero pairs.

% Find negligible residues:
ii = find(abs(res) < 1e-13 * norm(F, inf));
ni = length(ii);
if ( ni == 0 )
    % Nothing to do.
    return
elseif ( ni == 1 )
    fprintf('1 Froissart doublet.\n')
else
    fprintf('%d Froissart doublets.\n', ni)
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = [];
    f(jj) = [];
end

% Remove support points z from sample set:
for jj = 1:length(z)
    Z(Z == z(jj)) = [];
    F(Z == z(jj)) = [];
end
m = length(z);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(F, 0, M, M);
Sf = diag(f);
C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
A = SF*C - C*Sf;                    % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
w = V(:,m);

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w);
[pol, res, zer] = prz(r, z, f, w);

end % End of CLEANUP().


%% Automated choice of sample set

function [r, pol, res, zer, zj, fj, wj, errvec] = ...
    aaa_autoZ(F, dom, tol, mmax, cleanup_flag, mmax_flag)
%

% Flag if function has been resolved:
isResolved = 0;

% Main loop:
for n = 5:14
    % Sample points:
    % Next line enables us to do pretty well near poles
    Z = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), 1 + 2^n).';
    [r, pol, res, zer, zj, fj, wj, errvec] = aaa(F, Z, 'tol', tol, ...
        'mmax', mmax, 'cleanup', cleanup_flag);
    
    % Test if rational approximant is accurate:
    reltol = tol * norm(F(Z), inf);
    
    % On Z(n):
    err(1,1) = norm(F(Z) - r(Z), inf);
    
    Zrefined = linspace(dom(1)+1.37e-8*diff(dom), dom(2)-3.08e-9*diff(dom), ...
        round(1.5 * (1 + 2^(n+1)))).';
    err(2,1) = norm(F(Zrefined) - r(Zrefined), inf);
    
    if ( all(err < reltol) )
        % Final check that the function is resolved, inspired by sampleTest().
        % Pseudo random sample points in [-1, 1]:
        xeval = [-0.357998918959666; 0.036785641195074];
        % Scale to dom:
        xeval = (dom(2) - dom(1))/2 * xeval + (dom(2) + dom(1))/2;
        
        if ( norm(F(xeval) - r(xeval), inf) < reltol )
            isResolved = 1;
            break
        end
    end
end

if ( ( isResolved == 0 ) && ~mmax_flag )
    warning('CHEBFUN:aaa:notResolved', ...
        'Function not resolved using %d pts.', length(Z))
end

end % End of AAA_AUTOZ().
