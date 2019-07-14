function w = subsref(obj, varargin)
%SUBSREF   Evaluate a BARYFUN (calls feval).

w = feval(obj, varargin{1}.subs{1});

end
