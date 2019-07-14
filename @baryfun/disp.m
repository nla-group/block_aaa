function disp(obj)
%DISP   Display information about a BARYFUN.

[m,n] = size(obj.Ck{1});
d = length(obj.zk)-1;
fprintf('\tBARYFUN object of block size %d-by-%d and degree %d.\n', m, n, d);

end 
