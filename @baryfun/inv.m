function Rinv = inv(R)
%INV   BARYFUN inversion.
%
%   Rinv = inv(R) returns a BARYFUN corresponding to the inverse of R.
%   Only works for square BARYFUNs.

[m,n] = size(R.Ck{1});
if m~=n
    error('Inversion works only for square BARYFUNs');
end

Rinv = R;
Rinv.Ck = R.Dk;
Rinv.Dk = R.Ck;
