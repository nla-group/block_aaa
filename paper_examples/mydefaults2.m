%% Setting up defaults.

close all hidden
clear all

randn('state',50);
rand('state',50);

height = 2/3;
set(0,...
'defaultfigureposition',[180 100 800 800*height],...
'defaultaxeslinewidth',2,...
'defaultaxesfontsize',22,...
'defaultlinelinewidth',3,...
'defaultpatchlinewidth',1,...
'defaultlinemarkersize',14,...
'defaulttextinterpreter','latex');

%{
try
    chebfun();
catch
    try 
        addpath('../../../../Matlab/chebfun_v4.2.2889/chebfun')
        disp('Chebfun has been added to your Matlab path.')
    catch
        disp('Please download Chebfun and add it to your Matlab path.')
    end
end
%}

