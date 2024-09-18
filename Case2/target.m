function L=target(x,observed_data)
%   x is the input values of d-dimensional parameters
%  observed is the observations of target variable
% L is the log transfer of the target distribution
%%  1. write the x into certain path txt file
fid=fopen('/group_homes/lzu_public/home/u120220909911/Summer/Last/LE/input_step.txt','w');
for i=1:size(x,2)
    fprintf(fid,'%24.16e\n',x(i));
end
fclose(fid);

%% 2. run model
system('/group_homes/lzu_public/home/u120220909911/Summer/Last/LE/run.exe');


%% 3. caculate L
model_data=importdata('/group_homes/lzu_public/home/u120220909911/Summer/Last/LE/output_LE.txt');
T=size(model_data);
sum_error=sum((model_data-observed_data).^2);
sigmav=sum_error/T;
L=(-T/2*(log(2*pi*sigmav))-sum_error/(2*sigmav));


