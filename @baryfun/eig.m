function evs = eig(obj)
%EIG    Return the nonlinear eigenvalues (roots) of a BARYFUN.
  
evs = nonlinear_eig(obj.Ck,obj.zk);

end
