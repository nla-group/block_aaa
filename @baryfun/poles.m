function xi = poles(obj)
%POLES    Return the poles of a BARYFUN.

xi = nonlinear_eig(obj.Dk,obj.zk);

end
