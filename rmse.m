function err = rmse(pts, F1, F2)
%RMSE Root mean squared error of the difference of two 
%     block functions F1 and F2 sampled at pts.
%
% RMSE = sqrt( 1/npts * sum_{i=1}^npts || F1(pts(i))-F2(pts(i) ||_F^2 )

npts = length(pts);
if iscell(F1)
    FF1 = F1;
else
    FF1 = cell(npts,1);
    for i = 1:npts
        FF1{i} = F1(pts(i));
    end
end
if iscell(F2)
    FF2 = F2;
else
    FF2 = cell(npts,1);
    for i = 1:npts
        FF2{i} = F2(pts(i));
    end
end

err = 0;
for i = 1:npts
    err = err + norm(FF1{i} - FF2{i},'fro').^2;
end
err = sqrt(err/npts);
