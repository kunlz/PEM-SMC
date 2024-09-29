
% ------------------------------------------------- %
% Target function of simulation for candidate model %
% ------------------------------------------------- %
function [L] = target(para_est, obs)

% -----  Model parameter estimation for Lorenz  model ---------------------

% dimesion of target observation
[m,n] = size(obs);

% Call candidate model : take the Lorenz model as the test
[G_model] = Lorenzfk_mex(para_est, [-10 10 25], 0.01, 20);

L = 0;
sum_error = zeros(1,n);
sigmav = zeros(1,n);

for i = 1 : n

    % sum of error & sigma for multiple/single obs.
    sum_error(i) = sum((G_model(:,i)-obs(:,i)).^2);
    sigmav(i) = sum_error(i)/m;

    % log transfer of the likelihood function: L
    L = L+(-m/2*(log(2*pi*sigmav(i)))-sum_error(i)/(2*sigmav(i)));
end
% -----  Model parameter estimation for Lorenz  model ---------------------


%%%%%% Target 1  : multidiemsional normmal distribution with 2 modes at each dimension [upto 32 D is useful]
% d=size(x,2);
% mu1=-5*ones(1,d);
% sigma1=eye(d);
% mu2=5*ones(1,d);
% sigma2=eye(d);
% f=1/3*mvnpdf(x,mu1,sigma1)+2/3*mvnpdf(x,mu2,sigma2);
% L=log(f);

%%%%%%% Target 2 : one dimension two mode normal distribution
% p = 0.4;
% mu =[-2, 10];
% sd=[1, 2];
% f=  p*normpdf(x,mu(1),sd(1))+(1-p)*normpdf(x, mu(2),sd(2));
% L=log(f);

%%%%%%% Target 3 : 2 dimension multimode normal distribution
% x = para_est;
% d=2;
% Nm=15;
% mu6=[-9 2];
% sigma6=eye(d);
% 
% mu7=[-8 -4];
% sigma7=eye(d);
% 
% mu3=[-8 9];
% sigma3=eye(d);
% 
% mu1=[-5 3];
% sigma1=eye(d);
% 
% mu5=[-2 -5];
% sigma5=eye(d);
% 
% mu11=[-2 10];
% sigma11=eye(d);
% 
% mu4=1*ones(1,d);
% sigma4=eye(d);
% 
% 
% mu8=[2 -9];
% sigma8=eye(d);
% 
% mu9=[3 7];
% sigma9=eye(d);
% 
% mu2=[5 1];
% sigma2=eye(d);
% 
% 
% mu10=[6 -7];
% sigma10=eye(d);
% 
% mu12=[10 0];
% sigma12=eye(d);
% 
% mu13=[10 7];
% sigma13=eye(d);
% 
% mu14=[13 -5];
% sigma14=eye(d);
% 
% mu15=[13 14];
% sigma15=eye(d);
% f=1/Nm*(mvnpdf(x,mu1,sigma1)+mvnpdf(x,mu2,sigma2)+mvnpdf(x,mu3,sigma3)+mvnpdf(x,mu4,sigma4)+...
%     mvnpdf(x,mu5,sigma5)+mvnpdf(x,mu6,sigma6)+mvnpdf(x,mu7,sigma7)+mvnpdf(x,mu8,sigma8)+...
%     mvnpdf(x,mu9,sigma9)+mvnpdf(x,mu10,sigma10)+mvnpdf(x,mu11,sigma11)+mvnpdf(x,mu12,sigma12)+...
%     +mvnpdf(x,mu13,sigma13)+mvnpdf(x,mu14,sigma14)+mvnpdf(x,mu15,sigma15));
% 
% L=log(f);


% -----  Model parameter estimation AWBM model
% L is the log transfer of the likelihood function
% global obs
% T=size(obs,1);
% G_model=AWBM(x);
% sum_error=sum((G_model-obs).^2);
% sigmav=sum_error/T;
% L=(-T/2*(log(2*pi*sigmav))-sum_error/(2*sigmav));

%-------------------------------------------------------------------------

end

