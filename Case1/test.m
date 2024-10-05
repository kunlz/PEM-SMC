
clc
clear

Np = 50;                      % number of particles
S = 2000;                       % number of iterations
bound = [-15,-5;
         15,15]; 

tic
[para_iter] = PEM_sampler(Np, S, bound); tic


