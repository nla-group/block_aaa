function s = size(obj, dim)
%SIZE    Return the size of a BARYFUN.
  
s = size(obj.Ck{1});
if nargin == 2,  s = s(dim); end

end