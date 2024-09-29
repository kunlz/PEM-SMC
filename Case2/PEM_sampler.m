% -------------------------------------------------------------- %
% PEM-SMC algorithm developed to sample from the target function %
% -------------------------------------------------------------- %
% Authors & Code: Dr. Gaofeng Zhu, Lanzhou University            %
%                 Dr. Kun Zhang,   University of Hong Kong       %
% Date:           Nov. 30, 2021                                  %
% -------------------------------------------------------------- %
function [para_iter] = PEM_sampler(Np, S, bound, obs)
% --------------
% function input
% Np    -- number of particles
% S     -- number of iterations
% bound -- parameter boundary
% --------------

% crossover propabilty
pc = 0.6;

% dimesion of the target function
dim = size(bound, 2);

% preloca the matric used to save samplers
para_iter = nan(Np, dim, S);

% Generate initial population from the sampling space
w = 1/Np .* ones(Np,1);
para_old = bound(1, :) + (bound(2, :) - bound(1, :)) .* rand(Np, dim);

% Prior distribution
% Assuming that the prior distribution of parameters is Uniform
L_P0 = 0;
LB = bound(1, :); % lower boundary
UB = bound(2, :); % upper boundary
for i = 1 : length(LB)
    L_P0 = L_P0 - log(UB(i) - LB(i));
end

% Identify the sequence {beta_s}
betas = sequencesGen(S);

% The last five values are forced to be 1
betas(S+1 : S+5) = 1;

% SMC Iteration --------------
para_iter(:, :, 1) = para_old;

accept = zeros(S,1);
cal_weight = zeros(Np,1);
% alp_old = zeros(Np,1);

reverseStr = ' --> ';
fprintf('PEM-SMC :: completed');
for stage = 2 : S
    
    % Iteration for Np particle
    betas_new = betas(stage);
    betas_old = betas(stage - 1);

    for k = 1 : Np  % Parfor 

        % log of target distribution
        alp = betas_new * target(para_old(k, :), obs);

        % log of prior distribution
        prix = (1 - betas_new) * L_P0;

        alp_old = betas_old * target(para_old(k, :), obs);
        prix_old = (1 - betas_old) * L_P0;
        cal_weight(k, 1) = (alp + prix) - (alp_old + prix_old);
    end

    w = w .* exp(cal_weight);
    w = w ./ (sum(w));

    % Resampling -- Section
    % Method II
    ind = ResampSys(w, Np);
    para_new = para_old(ind, :);
    w = 1/Np .* ones(Np,1);

    % MCMC Step ---------------
    % rand walk
    b = 10^(-2).*ones(1,dim);

    for i = 1 : Np
        e = -b + 2.*b.*rand(1, dim);
        para_median = para_old(i, :) + e;
        alp_old = betas_new .* target(para_new(i, :), obs);
        alp_media = betas_new .* target(para_median, obs);      
 
        % M-H accept
        ratio = min([1, exp(alp_media - alp_old)]);
        u = rand;
        if u < ratio
            para_new(i, :) = para_median;  
            accept(stage) = accept(stage) + 1;
        else
            para_new(i, :) = para_new(i, :); 
        end        
    end    

    % Crossover Operater
    [para_median, index] = crossover(para_new, pc);

    % Iteration for Np particle
    for k = 1 : Np/2
        
        % Crossover step
        old_ind = index(2*k-1 : 2*k); 
        alp_old = betas_new * (target(para_new(old_ind(1), :), obs)...
            + target(para_new(old_ind(2), :), obs));

        alp_media = betas_new * (target(para_median(2*k-1, :), obs)...
            + target(para_median(2*k, :), obs));

        % M-H accept
        ratio = min([1, exp(alp_media - alp_old)]);
        u = rand;

        if u < ratio
            para_new(old_ind(1), :) = para_median(2*k-1, :);
            para_new(old_ind(2), :) = para_median(2*k, :);
            accept(stage) = accept(stage) + 1;
        else
            para_new(old_ind(1), :) = para_new(old_ind(1), :); % parfor
            para_new(old_ind(2), :) = para_new(old_ind(2), :);
        end

    end

    % MCMC Step ---------------
    % Mutation Operater
    for k = 1 : Np

        para_median2 = Generatep(para_new, k, bound(1, :), bound(2, :));
        alp_old = betas_new * target(para_new(k, :), obs);
        alp_media = betas_new * target(para_median2, obs);

        % M-H accept
        ratio = min([1, exp(alp_media - alp_old)]);
        u = rand;

        if u < ratio
            para_new(k, :) = para_median2;
            % accept(stage) = accept(stage) + 1;
        else
            para_new(k, :) = para_new(k, :);
        end

    end

    % Replace & Save Parameters
    para_old = para_new;
    para_iter(:, :, stage) = para_old;
    
    % ProgressBar
    percentDone = 100 * stage / S;
    msg = sprintf('%3.2f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]); 
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
disp(' --> Finish ')

end

% -------------------------------------------------------------------------
% Sub-functions -- 1
function [betas] = sequencesGen(S)
% --------
% function input:
% S - the number of generation
% --------

% Method I
% Expoential {beta_s} sequence (E. Jeremiah et al., 2012)
ss = 0 : S;
ss = ss';
m = 7e-8;

% using the relationship: phas=m*s^x to determine x
x = log(1/m) / log(S);
betas = m * ss.^x;

% Method II
% s=0:S;
% s=s';
% betas=s./S;

end

% -------------------------------------------------------------------------
% Sub-functions -- 2
function [new_gen, indexx] = crossover(old_para, pc)
% Crossover step
[~, indexx] = sort(rand(size(old_para, 1), 1));
match_gene = old_para(indexx, :);
pair = size(match_gene, 1) / 2;
% b = 1e-6;
% dem = size(match_gene, 2);
bits = size(match_gene, 2);
capirs = rand(pair, 1) < pc;
cpoints = randi(bits, [pair, 1]);
cpoints = cpoints .* capirs;

for i = 1 : pair
    new_gen([2*i - 1, 2*i], :) = [match_gene([2*i - 1, 2*i],...
        1:cpoints(i)), match_gene([2*i, 2*i - 1], cpoints(i) + 1:end)];
end

end

% -------------------------------------------------------------------------
% Sub-functions -- 3
function [para_new] = Generatep(para_old, k, LB, UB)
% Propose: generate a candidate partilce for MCMC mutation procedure using
%          the different evolution Markov Chain algorithm (ter Braak, 2006)
% --------

% Method I: Ter Braak, DE-MC
[N, d] = size(para_old);
gamma = 2.38 / sqrt(2*d);
b = 1e-6;

para_new = zeros(1, d);
while 1

    % select two particle that is different current particle k
    while true

        R1 = floor(rand * N);
        if R1 ~= k && R1 > 0

            break;
        end

    end

    while true

        R2 = floor(rand * N);
        if R2 ~= R1 && R2 ~= k && R2 > 0

            break
        end

    end

    % generate new canditate
    for j = 1 : d
        para_new(1, j) = para_old(k, j) +...
            gamma * (para_old(R1, j) - para_old(R2, j)) + b * randn;
    end

    lbjuest = para_new > LB;
    ubjuest = para_new < UB;

    if sum(lbjuest) == d && sum(ubjuest) == d

        break
    end

end

end

% -------------------------------------------------------------------------
% Sub-functions -- 4
function [outIndex] = ResampSys(w, N)
% Draw a total of N samples with probabilities proportional to the weight
% vector w, using Systematic Resampling algorithm.
%
% w            : normalized weight vector (sum to one)
% N (optional) : total number of samples; default to length(w)
% outIndex     : each element is an index into w, or, the "parent" of
%                the sample. Therefore if {X, w} is the original
%                particles-weights pair, then {X(outIndex), 1/N}
%                will be the resampled pair.
%
% from: Lingji Chen

eps = 1e-8; % small but not too small
len = length(w);
F = cumsum(w);

if abs(F(end) - 1) > eps
    error('the weight vector should be normalized.');
end

switch nargin

    case 1
        N = len;

    case 2
    otherwise
        error('wrong number of arguments');
end

s = rand / N;
inc = 1 / N;

outIndex = zeros(1, N);
j = 1;

for i = 1 : N

    while F(j) < s
        j = j + 1;
    end

    outIndex(i) = j;
    s = s + inc;
end

end
% ------------------------- E-N-D -----------------------------------------