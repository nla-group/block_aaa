%BARYFUN    BARYFUN constructor.
%
% Usage:
% obj = baryfun(zk, Ck, Dk) constructs a BARYFUN from the
%    support points in zk and the numerator/denominator 
%    coefficient matrices (cell arrays) in Ck/Dk. 

classdef baryfun
    properties
    zk
    Ck
    Dk
    transpose
    end % properties

    methods(Access = public, Static = true)
        % Constructor. 
        function obj = baryfun(varargin)

          obj.zk = varargin{1};
          obj.Ck = varargin{2};
          obj.Dk = varargin{3};
          if length(varargin) < 4
              obj.transpose = 0;
          else
              obj.transpose =  varargin{4};
          end

        end % function baryfun
    end % methods    
end % classdef baryfun
