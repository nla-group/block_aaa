function [Rbary,rmse,zk,Y] = baryfit5(F,lam,ell,opts)
%BARYFIT5   - fitting based on matrix-valued vector fitting (VF) in:
% https://www.sintef.no/projectweb/vectfit/downloads/matrix-fitting-toolbox/
% the rational matrix function Rbary is written in  

if nargin < 4
   opts.cmplx_ss=1;
   opts.poletype='linlogcmplx'; %Mix of linearly spaced and logarithmically spaced poles
   opts.Niter1=0;    %Number of iterations for fitting sum of elements (fast!) --> Improved initial poles
   opts.Niter2=3;    %Number of iterations for matrix fitting 
   opts.asymp=2;     %Fitting includes D but not E 
   opts.plot=0;
   opts.cmplx_ss=1;
   opts.stable = 0;
   opts.weightparam=1;
   opts.screen=0;
end

poles=[]; 
opts.N=ell ;%           %Order of approximation. 

p = length(lam);

% if F is function handle, evaluate it at all the lam's
if iscell(F)
    for i = 1:p
        FF(:,:,i) = F{i};
    end
else
    for i = 1:p
        FF(:,:,i) = F(lam(i));
    end
end
     
SER = VFdriver(FF,lam,poles,opts); %Creating state-space model and pole-residue model 

zk = SER.poles;

    % function handle to barycentric evaluation
    Rbary = @(z) SER.D+z*SER.E; 
    for j = 1:ell
        Rbary =  @(z) Rbary(z) + SER.R(:,:,j)/(z-SER.poles(j));
    end
    
 
%  
%  
% % the set of residues Y in cell format
% Y = cell(ell,1);
% for jj = 1:ell    
%     Y{jj} = SER.R(:,:,j);
% end
%     
err = zeros(1,p);

% compute rmse
for i = 1:p
   err(i) = norm(FF(:,:,i)-Rbary(lam(i)),'fro').^2; 
end
rmse = sqrt(sum(err)/p);


end

