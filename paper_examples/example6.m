% compare two fitting algorithms (RKFIT5 and Block-AAA) for a simple
% low-order (rational) function
%
% with addition of low, moderate, and hig levels of normally distributed NOISE 

addpath('matrix_fitting_toolbox_1')
addpath('block_aaa')
addpath('rktoolbox')
clc
mydefaults2
warning('off','all')

npts = 500;  % number of sampling points
pts = logspace(-1,1,npts)*1i;

%% option 1
F = @(z) (z-1)/(z^2+z+2);

% %% option 2
% F = @(z) [ (z^2 - 3)/(z^2+6*z+5) 1/(z+10) ;
%    z/(z^3+3*z+1) z/10 ];
% 
% %% option 3 (20th order Chebyshev type I filter)
% k = 10;
% [A,B,C,D] = cheby1(k,3,[500 560]/750);
% C = 10^5*C;
% B = 10^5*B;
% F = @(z) C*((z*eye(2*k)-A)\B)+D;

nlevel = 10^(-2);

FF = {}; 
noise = [];
for i = 1:npts
    noise(i) = nlevel*randn(1)+nlevel*randn(1)*1i;
    FF{i} = F(pts(i))+noise(i);
end

[m,n] = size(FF{1});

maxell = 5; % maximum degree tried for all methods
timingell = 5; % degree for which to time each of the methods
runs = 1; % number of runs for the timings

%% plot functions
figure(1)
% m=n=1 here
   for i = 1:m
      for j = 1:n
         for l = 1:npts 
             Fval(l) = FF{l}(i,j);
         end
         loglog(imag(pts),abs(Fval),'k-'); hold on;
      end
   end
    

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
    ops.iter = 3;
    ops.tol = 0;
    xi = inf(1,ell);
    rng('default')
    [R,rmse_] = baryfit6(FF,pts,xi,ops);
    ERR4(ell+1) = rmse_;
end
toc

figure(1)
hold on
Rval = arrayfun(@(x) R(x),pts);
loglog(imag(pts),abs(Rval),'r-');

figure(2)
loglog(imag(pts),abs(Fval - Rval)+eps,'r-')
xlabel('$z/i$')
ylabel('error')
ylim([1e-4,1])
shg

%%
ops.iter = 3;
ops.tol = 0;
xi = inf(1,timingell);
tic
for r = 1:runs
    rng('default')
    [~,rmse_] = baryfit6(FF,pts,xi,ops);
end
TIME4 = 1000*toc/runs;

% % fixed number of RKFIT iterations with full initial poles at inf
% disp('standard rkfit (fixed iters) -- note that we are doing superfluous work!')
% tic
% for ell = 0:maxell
%     disp(['ell=' num2str(ell)])
%     ops.iter = 10;
%     ops.tol = 0;
%     xi = inf(1,ell);
%     rng('default')
%     [R,rmse_] = baryfit6(FF,pts,xi,ops);
%     ERR44(ell+1) = rmse_;
% end
% toc
% 
% ops.iter = 10;
% ops.tol = 0;
% xi = inf(1,timingell);
% tic
% for r = 1:runs
%     rng('default')
%     [R,rmse_] = baryfit6(FF,pts,xi,ops);
% end
% TIME44 = 1000*toc/runs;

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

%%
figure(1)
hold on
Rval_aaa = arrayfun(@(x) Rbary(x),pts);
loglog(imag(pts),abs(Rval_aaa),'b-');
legend('function','RKFIT','AAA')
xlabel('$z/i$')
ylabel('$|f(z)|$')
title(['scalar example with noisy data'])
shg

figure(2)
hold on
loglog(imag(pts),abs(Fval - Rval_aaa)+eps,'b-')
loglog(imag(pts),abs(noise)+eps,'k:','Color',[.5,.5,.5])
legend('RKFIT','AAA','true noise','Location','SouthEast')
xlabel('$z/i$')
ylabel('$|f(z) - r(z)|$')
title(['scalar example with noisy data'])
shg
ylim([1e-5,1])
set(gca,'YTick',[1e-4,1e-2,1])


%%
warning('off','all')
ops.maxit = timingell+1;
ops.tol = 0;
tic
for r = 1:runs
    rng('default')
    [R,rmse_] = block_aaa(FF,pts,ops);
end
TIME5 = 1000*toc/runs;

% %%
% % plot poles and roots of approximant
% %[Rbary,ERR2_,zk,Ck,Dk] = baryfit2(FF,pts,5);
% %poles = nonlinear_eig(Dk,zk);
% %rts = nonlinear_eig(Ck,zk);
% %figure; plot(pts,'k.'); hold on; plot(poles,'ro'); plot(rts,'b+')
% 
% % noninterpolatory block-FIT
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


%%
clc
disp(['ERRORS AND TIMINGS FOR DEGREE ' num2str(timingell) ' APPROXIMANTS'])

figure(3)

disp('RKFIT 5 iter')
semilogy(0:maxell,ERR4,'r-o'); % RKFIT
hold on
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR4(timingell+1),TIME4))

disp('block-AAA')
semilogy(0:maxell,ERR5,'b-o'); % block-AAA
disp(sprintf('    error %10.3d   --    timing %5.1f ms',ERR5(timingell+1),TIME5))
% 
% yt={'10^{-5}'; '10^{-4}';'10^{-3}'; '10^{-2}'; '10^{-1}';'10^{0}'}; 
% set(gca,'ytick',[10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^(0)]); 
% set(gca,'yticklabel',yt);
title(['scalar example with noisy data'])
axis tight
ylim([1e-2,2])
xlabel('degree d')
ylabel('RMSE')
grid on
legend('RKFIT','AAA','Location','NorthEast');
