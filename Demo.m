%% Test 1
% ------------------------------------------- E・ Lorenz model
clear
clc

load dataLorenz % target/observed data

obs = dataLorenz;
Np = 100;                      % number of particles
S = 700;                       % number of iterations

% Lorenz model parameter calibration 
% Parameter boundary
%       deta belta pho
bound = [6,  1,  20;
         14, 6,  36]; 

tic
[paramter_iteration] = PEM_sampler(Np, S, bound, obs); toc

% tic
% [paramter_iteration] = PEM_sampler_o(Np,S,bound,obs);toc

%% 
clear paramete_trajectory
for i = 1 : S
    parameter = paramter_iteration(:,:,end);  
    paramete_trajectory(i,:) = paramter_iteration(12,:,i);
end

% parameter used to generate observation
para = [10 8/3 28]; 

% plot(para,median(parameter),'.')
for i=1:length(para)
    plot(paramete_trajectory(:,i))
    hold on
    plot(S,para(i),'x')
end
para_est = median(parameter);

%% figure(1)
close all
[x3] = Lorenzfk_mex(para_est,[-10 10 25],0.01,20);

plot(x3(:,2))
hold on
plot(obs(:,2))
%%
close all
plot(x3(:,1))
hold on
plot(obs(:,1))
%%
close all
plot(obs(:,1),obs(:,3),'LineWidth',3)
hold on
plot(x3(:,1),x3(:,3),'LineWidth',1.5)
