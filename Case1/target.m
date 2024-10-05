function L = target(x)
% x is the sample
% L is the log transfer of the target distribution

%--------------------------------------------------------------------------
d=2;
Nm=15;


mu1=[-5 3];
sigma1=eye(d);

mu2=[5 1];
sigma2=eye(d);

mu3=[-8 9];
sigma3=eye(d);

mu4=1*ones(1,d);
sigma4=eye(d);

mu5=[-2 -5];
sigma5=eye(d);

mu6=[-9 2];
sigma6=eye(d);

mu7=[-8 -4];
sigma7=eye(d);

mu8=[2 -9];
sigma8=eye(d);

mu9=[3 7];
sigma9=eye(d);

mu10=[6 -7];
sigma10=eye(d);

mu11=[-2 10];
sigma11=eye(d);

mu12=[10 0];
sigma12=eye(d);

mu13=[10 7];
sigma13=eye(d);

mu14=[13 -5];
sigma14=eye(d);

mu15=[13 14];
sigma15=eye(d);


f=1/Nm*(mvnpdf(x,mu1,sigma1)+mvnpdf(x,mu2,sigma2)+mvnpdf(x,mu3,sigma3)+mvnpdf(x,mu4,sigma4)+...
    mvnpdf(x,mu5,sigma5)+mvnpdf(x,mu6,sigma6)+mvnpdf(x,mu7,sigma7)+mvnpdf(x,mu8,sigma8)+...
    mvnpdf(x,mu9,sigma9)+mvnpdf(x,mu10,sigma10)+mvnpdf(x,mu11,sigma11)+mvnpdf(x,mu12,sigma12)+...
    +mvnpdf(x,mu13,sigma13)+mvnpdf(x,mu14,sigma14)+mvnpdf(x,mu15,sigma15));

L=log(f);

end
