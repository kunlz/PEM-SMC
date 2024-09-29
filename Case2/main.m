clc
clear all
%% PEM-SMC������

Np=50;                      % number of particles
S=300;                      % number of iterations
% %  Parameter boundary
% %   vmax25  ��    m    ��_s    b      ��_(sat-d)
bound=[10,  0.06,  4,   50,    0.01,   0.35;
       200, 0.08,  9,   500,   0.04,   0.55];  
[paramter_iteration]=PEM_sampler(Np,S,bound); 
save parameter paramter_iteration