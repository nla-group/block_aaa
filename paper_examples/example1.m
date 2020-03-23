addpath('block_aaa')
addpath('rktoolbox')
clc
mydefaults2
warning('off','all')

% % first choice of an academical example
% F = @(z) [ (z^2 - 3)/(z+1) , 1 , 1/(z-4) ;
%     z/(z^3+1) , z^2 , 7 ;
%     1/(z^3 + z^2) , z^2-8 , z ];

% determines whether F is symmetric or not
sym = 0;

if(sym)
    F = @(z) [ 2/(z+1) (3-z)/(z^2-5+z) ;
        (3-z)/(z^2-5+z) (2+z^2)/(z^3 + 3*z^2-1)];
else
    F = @(z) [ 2/(z+1) (3-z)/(z^2+5+z) ;
        (3-z)/(z^2-5+z) (2+z^2)/(z^3 + 3*z^2-1)];
end

npts = 100;
%pts = 2*exp(2i*pi*(1:npts)/npts);        % sampling points
pts = logspace(0,2,npts)*1i;

FF = {}; 
for i = 1:npts
    FF{i} = F(pts(i));
end

[m,n] = size(FF{1});


maxell = 20; % maximum degree tried for all methods
timingell = 10; % degree for which to time each of the methods
runs = 20; % number of runs for the timings

%% plot functions
figure

   for i = 1:m
      for j = 1:n
         for l = 1:npts 
             Fval(l) = FF{l}(i,j);
         end
         loglog(imag(pts),abs(Fval)); hold on;
      end
   end
   
xlabel('$z/i$')
ylabel('$|F_{i,j}(z)|$')
if(sym)
   title('symmetric toy example')
else
   title('nonsymmetric toy example')
end


%% set-valued AAA (this is inefficient; TODO: change aaa_sv to return RMSE)
disp('set-valued AAA')
tic
for ell = 0:maxell
    rng('default')
    [R,rmse_] = baryfit3(FF,pts,ell);
    ERR1(ell+1) = rmse_;
end
toc

tic
for r = 1:runs
    rng('default')
    [R,rmse_] = baryfit3(FF,pts,timingell);
end
TIME1 = 1000*toc/runs;

%% surrogate AAA (this is inefficient; TODO: change aaa to return RMSE)
disp('surrogate AAA')
tic
for ell = 0:maxell
    rng('default')
    [R,rmse_] = baryfit4(FF,pts,ell);
    ERR2(ell+1) = rmse_;
end
toc

tic
for r = 1:runs
    rng('default')
    [R,rmse_] = baryfit4(FF,pts,timingell);
end
TIME2 = 1000*toc/runs;

%% vector fitting 5 iterations
opts.Niter2=5; 
opts.cmplx_ss=1; opts.poletype='linlogcmplx'; opts.Niter1=0; opts.asymp=2;   
opts.plot=0; opts.cmplx_ss=1; opts.stable = 0; opts.weightparam=1; opts.screen=0;

disp('vector fitting')
vfit_tic = tic;
for ell = 0:maxell
    if ell == 0
        ERR3(1) = NaN;
        continue
    end
    if(norm(imag(pts)*1i-pts)&&mod(ell,2))
        rmse_ = 1;
    else 
        rng('default')
        [R,rmse_] = baryfit5(FF,pts,ell,opts);
    end
    ERR3(ell+1) = rmse_;
end
toc(vfit_tic)

warning('off','all')
vfit_tic = tic;
for r = 1:runs
    rng('default')
    [R,rmse_] = baryfit5(FF,pts,timingell,opts);
end
TIME3 = 1000*toc(vfit_tic)/runs;


%% vector fitting 10 iterations
opts.Niter2=10; 

disp('vector fitting')
vfit_tic = tic;
for ell = 0:maxell
    if ell == 0
        ERR33(1) = NaN;
        continue
    end
    if(norm(imag(pts)*1i-pts)&&mod(ell,2))
        rmse_ = 1;
    else 
        rng('default')
        [R,rmse_] = baryfit5(FF,pts,ell,opts);
    end
    ERR33(ell+1) = rmse_;
end
toc(vfit_tic)

warning('off','all')
vfit_tic = tic;
for r = 1:runs
    rng('default')
    [R,rmse_] = baryfit5(FF,pts,timingell,opts);
end
TIME33 = 1000*toc(vfit_tic)/runs;

%% incremental rkfit
% disp('incremental rkfit - note that we are doing superfluous work')
% tic
% for ell = 0:maxell 
%     disp(['ell=' num2str(ell)])
%     ops.iter = ell;
%     ops.tol = 0;
%     [R,rmse_] = baryfit6(FF,pts,[],ops);
%     ERR4_inc(ell+1) = rmse_;
% end
% toc

%% fixed number of RKFIT iterations with full initial poles at inf
disp('standard rkfit (fixed iters) -- note that we are doing superfluous work!')
tic
for ell = 0:maxell
    %disp(['ell=' num2str(ell)])
    ops.iter = 5;
    ops.tol = 0;
    xi = inf(1,ell);
    rng('default')
    [R,rmse_] = baryfit6(FF,pts,xi,ops);
    ERR4(ell+1) = rmse_;
end
toc

ops.iter = 5;
ops.tol = 0;
xi = inf(1,timingell);
tic
for r = 1:runs
    rng('default')
    [R,rmse_] = baryfit6(FF,pts,xi,ops);
end
TIME4 = 1000*toc/runs;

%% fixed number of RKFIT iterations with full initial poles at inf
disp('standard rkfit (fixed iters) -- note that we are doing superfluous work!')
tic
for ell = 0:maxell
    %disp(['ell=' num2str(ell)])
    ops.iter = 10;
    ops.tol = 0;
    xi = inf(1,ell);
    rng('default')
    [R,rmse_] = baryfit6(FF,pts,xi,ops);
    ERR44(ell+1) = rmse_;
end
toc

ops.iter = 10;
ops.tol = 0;
xi = inf(1,timingell);
tic
for r = 1:runs
    rng('default')
    [R,rmse_] = baryfit6(FF,pts,xi,ops);
end
TIME44 = 1000*toc/runs;

%% interpolatory block-AAA
disp('interpolatory block-AAA')
% tic
% [Rbary,ERR5,zk,Ck,Dk] = baryfit2(FF,pts,maxell);
% toc
warning('off','all')
opts.tol = 0;
opts.maxit = maxell+1;
opts.return = 'last';
tic
rng('default')
[Rbary,ERR5,out1] = block_aaa(FF,pts,opts); 
toc

warning('off','all')
ops.maxit = timingell+1;
ops.tol = 0;
tic
for r = 1:runs
    rng('default')
    [R,rmse_] = block_aaa(FF,pts,ops);
end
TIME5 = 1000*toc/runs;

%%
% plot poles and roots of approximant
%[Rbary,ERR2_,zk,Ck,Dk] = baryfit2(FF,pts,5);
%poles = nonlinear_eig(Dk,zk);
%rts = nonlinear_eig(Ck,zk);
%figure; plot(pts,'k.'); hold on; plot(poles,'ro'); plot(rts,'b+')

% %% noninterpolatory block-FIT
% opts.svd = 1;
% opts.scale = 0;
% opts.plot = 0;
% opts.iter = 1;
% disp('noninterpolatory block-FIT')
% tic
% for ell = 0:maxell
%     rng('default')
%     [Rbary,err] = baryfit1(FF,pts,ell,opts);
%     ERR6(ell+1) = max(err);
% end
% toc
% 
% tic
% for r = 1:runs
%     rng('default')
%     [R,rmse_] = baryfit1(FF,pts,timingell,opts);
% end
% TIME6 = 1000*toc/runs;

%% Loewner 
for ell = 0:maxell
    rng('default')
    [R,rmse_] = loewfit3(FF,pts,ell);
    ERR7(ell+1) = rmse_;
end

tic
for r = 1:runs
    rng('default')
    [R,rmse_] = loewfit3(FF,pts,timingell);
end
TIME7 = 1000*toc/runs;

%%
clc
disp(['ERRORS AND TIMINGS FOR DEGREE ' num2str(timingell) ' APPROXIMANTS'])
figure
disp('set-valued AAA')
semilogy(0:maxell,ERR1,'b-'); % set-valued AAA
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR1(timingell+1),TIME1))
hold on
xlabel('order d')
ylabel('RMSE')
grid on
axis([0,maxell,1e-16,10])
%axis tight
yt={'10^{-15}'; '10^{-10}'; '10^{-5}'; '10^{0}'} ; 
set(gca,'ytick',[10^(-15) 10^(-10) 10^(-5) 10^(0)]); 

disp('surrogate AAA')
semilogy(0:maxell,ERR2,'c:','Color',[.3,.5,1]); % surrogate AAA
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR2(timingell+1),TIME2))

disp('VF 5 iter')
semilogy(0:maxell,ERR3,'g-','Color',[0,0.7,0]); % VF 
hold on;
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR3(timingell+1),TIME3))

disp('VF 10 iter')
semilogy(0:maxell,ERR33,'m--'); % VF
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR33(timingell+1),TIME33))


disp('RKFIT 5 iter')
semilogy(0:maxell,ERR4,'g','Color',[0.59,0.29,0]); % RKFIT
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR4(timingell+1),TIME4))

disp('RKFIT 10 iter')
semilogy(0:maxell,ERR44,'c--'); % RKFIT
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR44(timingell+1),TIME44))


disp('Loewner')
semilogy(0:maxell,ERR7,'g:','Color',[1,0.6,0.6]); % Loewner
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR7(timingell+1),TIME7))

disp('block-AAA')
semilogy(0:maxell,ERR5,'r-'); % block-AAA
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR5(timingell+1),TIME5))


% disp('block-NONFIT')
% h6=semilogy(0:maxell,ERR6,'r:','Color',[1,.6,.6]); % block-NONINT
% disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR6(timingell+1),TIME6))

legend('set-valued AAA','surrogate AAA','VF (5 iter)','VF (10 iter)','RKFIT (5 iter)','RKFIT (10 iter)','Loewner','block-AAA','Location','NorthEast');

if(sym)
   title('symmetric toy example')
else
   title('nonsymmetric toy example')
end
