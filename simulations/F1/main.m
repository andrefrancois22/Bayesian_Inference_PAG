% ==> Ideal bayesian Observer Simulation
clear all; close all; clc;
% ==> lib
addpath('pfc_functions/');

rng(1,"twister");

% ==> params
% ======> number of orientations
n_o = 7; 
% ==> number of orientations must be odd for this simulation! 
% e.g. include an ambiguous vertical stimulus (identical to experiment)
assert(mod(n_o,2),'>>>>>> the number of orientations must be odd! <<<<<<')

% ======> sigmas (sensory noise range)
sigs = linspace(0.25, 3.75, 7); 

% ======> (1) define a prior pr probability ratio
fx = 2; % (experiments have a 2x factor)

% => central (vertical) orientation index (don't change this)
vidx = round(n_o/2);
if ~mod(n_o,2)
    vidx = vidx + 1;
end

% => option summing truncated likelihood density and adding to extreme orientations
% => 'true'  means keep truncated (normalized) densities 
% => 'false' means add the truncated Gaussian densities to edges 
% LK_TRC_FLAG = 'true';
LK_TRC_FLAG = 'false';

% ======> number of 'sample' trials per orientation (Robbe's suggestion)
% => lower case k is just the number of repeat samples for each orientation
% => grid (rather than random sample)
k = 1000; %1000 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> form the prior (step function)
pr = ones(1,n_o);
pr(vidx:end) = pr(vidx:end)*fx;
% => normalization & check
pr = pr./sum(pr); assert(sum(pr));

% ==> two contexts
ctxs = {'cw','ccw'};

% For storing accuracy values
mle_v = nan(length(ctxs),length(sigs));
map_v = nan(length(ctxs),length(sigs));

% For storing proportions and counts by orientation 
pr_mle = nan(length(ctxs),2,n_o); 
pr_map = nan(length(ctxs),2,n_o);

% ==> for storing counts (for PF model fit)
% => cw
cw_ct_mle = nan(length(ctxs),length(sigs),n_o);  
cw_ct_map = nan(length(ctxs),length(sigs),n_o);
% => ccw
ccw_ct_mle = nan(length(ctxs),length(sigs),n_o); 
ccw_ct_map = nan(length(ctxs),length(sigs),n_o);

% ==> for each context
for I = 1:2
    
    % => context    
    ctx = ctxs{I};
    
    % ==> if ccw context, switch prior
    if strcmp('cw',ctx)
        s0 = [repmat([1:(vidx-1)]',1*k,1); repmat([vidx:n_o]',fx*k,1)]';
        % ==> (2) total number of K stimulus orientation samples (from prior)
        K = length(s0);          
    elseif strcmp('ccw',ctx)
        pr = fliplr(pr);
        s0 = [repmat([1:(vidx)]',fx*k,1); repmat([(vidx+1):n_o]',k,1)]';
        % ==> (2) total number of K stimulus orientation samples (from prior)
        K = length(s0);          
    end  
    
    i = 1;
    for sig = sigs
        % ==> compute accuracy for map, mle, density estimates etc
        % => return liklihoods and posteriors
        [mle,map,pos,acc_mle,acc_map,lk,po] = simpredf(pr,s0,vidx,sig,n_o,K,LK_TRC_FLAG);    
        % ==> store accuracy values
        mle_v(I,i) = acc_mle; map_v(I,i) = acc_map; 

        % ==> for PF curves and proportions & counts
        for or = 1:n_o
            
            % ==> compute counts 'cw' (e.g. strictly greater than vertical)
            cw_ct_mle(I,i,or)  = sum( (mle(s0' == or) > vidx) & s0(s0' == or)' );
            cw_ct_map(I,i,or)  = sum( (map(s0' == or) > vidx) & s0(s0' == or)' );
            % ==> compute counts 'ccw' (e.g. strictly less than vertical)
            ccw_ct_mle(I,i,or) = sum( (mle(s0' == or) < vidx) & s0(s0' == or)' );
            ccw_ct_map(I,i,or) = sum( (map(s0' == or) < vidx) & s0(s0' == or)' );          

            % ==> compute proportion 'cw' (e.g. strictly greater than vertical)            
            pr_mle(I,i,or) =  cw_ct_mle(I,i,or) / sum(s0' == or);
            pr_map(I,i,or) =  cw_ct_map(I,i,or) / sum(s0' == or);        
        end
        % index for sigma 
        i = i + 1;    
    end
end

fprintf('computing psychometric functions...\n')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> PF curve estimation Set options
fitFlag         = 1;       % 0 or 1, where 1 means do the fit right now
nBootStraps     = 100;     % determines number of non parametric bootstraps if fitFlag == 1
options         = optimoptions('fmincon');
options.Display = 'off';
% Set bounds on model parameters
% Model fit to choice data split by task context (1 & 2: guess rate; 3: perceptual uncertainty; 4 & 5: decision criterion)
startVec_M1 = [0.04 0.04 1 0 0];
startVec_V2 = [0.04 0.04 3 -15 15];
LB_M1(1,1)  = 0;                         UB_M1(1,1) = 0.05;      % lapse rate
LB_M1(2,1)  = 0;                         UB_M1(2,1) = 0.05;      % lapse rate
LB_M1(3,1)  = 0.1;                       UB_M1(3,1) = 10;         % perceptual uncertainty
LB_M1(4,1)  = -200;                      UB_M1(4,1) = 200;        % decision criterion (-10->10 range originally)
LB_M1(5,1)  = -200;                      UB_M1(5,1) = 200;        % decision criterion
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% orientations - for 7 this becomes { -3, -2, -1, 0, +1, +2, +3 }
d = (1:n_o) - vidx;

% => bias vectors
bsv_map = nan(length(sigs),1);
bsv_mle = nan(length(sigs),1);

for i = 1:length(sigs)
    
    % ==> counts for MAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> context 1 ccw response (incongruent) & ==> context 2 ccw response (congruent)
    nccw_map = [squeeze(ccw_ct_map(1,i,:))';  squeeze(ccw_ct_map(2,i,:))']; 
    % ==> context 1 cw response (congruent)    & ==> context 2 cw response (incongruent)
    ncw_map  = [squeeze(cw_ct_map(1,i,:))';   squeeze(cw_ct_map(2,i,:))' ]; 
    % ==> pass through PF function fit
    obFun    = @(paramVec) giveNLL(paramVec, d, nccw_map, ncw_map, 'M1');         
    
    if i < .75*numel(sigs)
        paramEst_M1_map = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);
    else
        paramEst_M1_map = fmincon(obFun, startVec_V2, [], [], [], [], LB_M1, UB_M1, [], options);
    end

    % ==> return pf
    [~, predPF_map] = giveNLL(paramEst_M1_map, d, nccw_map, ncw_map, 'M1');
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    
    % ==> counts for MLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ==> context 1 ccw response (incongruent) & ==> context 2 ccw response (congruent)
    nccw_mle = [squeeze(ccw_ct_mle(1,i,:))';  squeeze(ccw_ct_mle(2,i,:))']; 
    % ==> context 1 cw response (congruent)    & ==> context 2 cw response (incongruent)
    ncw_mle  = [squeeze(cw_ct_mle(1,i,:))';   squeeze(cw_ct_mle(2,i,:))' ];     
    % ==> pass through PF function fit
    obFun    = @(paramVec) giveNLL(paramVec, d, nccw_mle, ncw_mle, 'M1');         
    paramEst_M1_mle = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);
    % ==> return pf
    [~, predPF_mle] = giveNLL(paramEst_M1_mle, d, nccw_mle, ncw_mle, 'M1');
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
    % ==> store bias estimates (map) - calculating as:
    % the distance between the points of subjective equality 
    % for the cw and ccw prior context PF curves
    bsv_map(i) = paramEst_M1_map(end) - paramEst_M1_map(end-1);   
    % store bias estimates (mle)
    bsv_mle(i) = paramEst_M1_mle(end) - paramEst_M1_mle(end-1);     
end


% ==> plot
for I = 1:2
    % => context    
    ctx = ctxs{I};
    % ==> plots
    pltf();
end
