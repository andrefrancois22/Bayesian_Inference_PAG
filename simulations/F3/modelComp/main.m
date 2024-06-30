clear all; close all; clc;
% ==> data directory
drc = '../../../data/';

% ==> dynamic range analysis types
dyn_types = {'popu','stim'};
% ==> load dynamic range results per stimulus orientation? or per population?
dyn_type = dyn_types{1};
% ==> load choice proportions and counts from dynamic range split on real data
% CPs (<Session ID>_<contrast>_<dynamic range split>)
load([drc,'CPs_dynr_',dyn_type,'.mat'],'CPs');
% CTs (<Session ID>_<contrast>_<dynamic range split>_<counts (nccw (1) and ncw (2))>)
load([drc,'CTs_dynr_',dyn_type,'.mat'],'CTs');
% ==> load PFs
% PFs (<Session ID>_<contrast>_<dynamic range split>)
load([drc,'PFs_dynr_',dyn_type,'.mat'],'PF');

% ==> parameters 
% => number of simulated trials
N = 500;
tm = 200; % timepoints (same resolution as actual DV fits)  
% ==> bound
bnd = 12; 
% ==> sd (for diffusion component)
sd = 1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> PF curve estimation Set options
fitFlag         = 1;       % 0 or 1, where 1 means do the fit right now
nBootStraps     = 10;     % determines number of non parametric bootstraps if fitFlag == 1
options         = optimoptions('fmincon','ConstraintTolerance', 1e-100, ...
                                         'OptimalityTolerance', 1e-100, ...
                                         'StepTolerance',       1e-100, ...
                                         'Algorithm', 'Interior-point');
options.Display = 'iter';
% Set bounds on model parameters
LB(1,1)  = 0.01;                      UB(1,1) = 0.06;    % orientation mean interval
LB(2,1)  = -0.2;                      UB(2,1) = 0.2;  % overall bias (not decision bias - this is left or right shift of all curves)
LB(3,1)  = 0;                         UB(3,1) = 20;   % drift rate (slope of the linear drift)
LB(4,1)  = 0;                         UB(4,1) = 200;   % initial offset
% LB(5,1)  = 6;                         UB(5,1) = 20;   % bound
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> for each session, and for each contrast level, fit acc-bound predicted dynr PFs curves to actual data
iS = 17;
cr = 2; % => lo contrast (fit models separately for different stimulus contrast images)
% ==> accumulation to bound model type
mod_type = 'm2';

% ==> fmincon optim often converges to bad local minima - sensitive to
% initial parameter settings. Initialize with best results from a coarse
% grid search

% ==> grid resolution - step
stp = 10;

itsvs = linspace(LB(1,1),UB(1,1),stp); % --> orientation mus interval
sftvs = linspace(LB(2,1),UB(2,1),stp); % --> horizontal shift
dftvs = linspace(LB(3,1),UB(3,1),stp); % --> drift rate
oftvs = linspace(LB(4,1),UB(4,1),stp); % --> initial offset

% ==> forward model NNLs
NNLs_init = nan(stp,stp,stp,stp);

% ==> brute force grid search for initial param settings before fmincon
i = 1;
for itvi = 1:stp                                    % --> orientation mus interval
    % => interval
    itv = itsvs(itvi);
    for sfti = 1:stp                                % --> horizontal shift
        % => shift
        sft = sftvs(sfti);
        for dfti = 1:stp                            % --> drift rate
            % => drift rate
            dft = dftvs(dfti);
            for ofti = 1:stp                        % --> initial offset
                % => offset
                oft = oftvs(ofti);               
                % ==> interval, overall bias (shift), fc (drift rate), initial offset
                startVec = [itv, sft, dft, oft];                
                % ==> run forward model    
                [nll, ~, ~] = modfitf(iS, cr, CTs, N, tm, sd, startVec, mod_type);
                % ==> store
                NNLs_init(itvi,sfti,dfti,ofti) = nll;
                % => for command line updates
                i = i + 1;
            end            
        end
    end
    % => update               
    clc; fprintf('completed step %d of %d...\n',i,stp^length(LB));  
end

% ==> return 4 indices of min NLL in NNLs_init
[min_nll, idx] = min(NNLs_init(:));
[i,j,k,l] = ind2sub( size(NNLs_init), idx );

% ==> interval, overall bias (shift), fc (drift rate), initial offset, bound
startVec = [itsvs(i),sftvs(j),dftvs(k),oftvs(l)]; %[0.038, 0, 8, 90];%, 12]; 

% ==> repeats
r = 100;

% ==> NLLs
nlls = cell(r,1);
% ==> simulation optim params
startVecs = cell(r,1);

% ==> objective function for fmincon
obFun = @(paramVec) modfitf(iS, cr, CTs, N, tm, sd, paramVec, mod_type); 

% ==> run fmincon several times (variable results even with multistart)
for rep = 1:r
    % ==> run optimization
    simfit = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
    % ==> return best fit simulated dynamic range proportions
    [nll, ~, ~] = modfitf(iS, cr, CTs, N, tm, sd, simfit, mod_type);
    % ==> track nlls
    nlls{rep} = nll;    
    % ==> reset startVec if new minimum for nll in next round
    % ==> kind of a simulated annealing procedure...
    if nll <= min(cell2mat(nlls))
        % ==> startVec
        startVec = simfit;
        % ==> store
        startVecs{rep} = simfit;
    end    
    fprintf('Finished simulated annealing step %d of %d...\n',rep,r);
end

% ==> find index of min
[~,I] = min(cell2mat(nlls));
% ==> set optimal paramters
simfit = startVecs{I};
% ==> run forward model
[nll, propcw, propccw] = modfitf(iS, cr, CTs, N, tm, sd, simfit, mod_type);
% ==> update
fprintf('Done with artificial multistart...\n')

% ==> cashe the results
save(['params/','iS_',num2str(iS),'_cr_',num2str(cr),'_startVec.mat'],'simfit','propcw','propccw','nll');

% ==> load optimal parameters and simulated choice proportions
ps = load(['params/','iS_',num2str(iS),'_cr_',num2str(cr),'_startVec.mat']);
propcw  = ps.propcw;
propccw = ps.propccw;
simfit  = ps.simfit; 

% ==> visualize the results
close all; 
figure(1);
subplot(2,3,1);
hold on; hold all;
plot(1:7,propcw(1,:),'r-h')
plot(1:7,propccw(1,:),'b-h')
plot(1:7,propcw(2,:),'rh:')
plot(1:7,propccw(2,:),'bh:')
subplot(2,3,2);
hold on; hold all;
plot(CPs{iS,cr,2}(2,:),'b:')
plot(CPs{iS,cr,2}(1,:),'r:')
plot(CPs{iS,cr,1}(1,:),'r')
plot(CPs{iS,cr,1}(2,:),'b')
subplot(2,3,3);
hold on; hold all;
plot(1:7, CTs{iS,cr,1,2}(2,:) ./ (CTs{iS,cr,1,2}(2,:) + CTs{iS,cr,1,1}(2,:)),'b');
plot(1:7, CTs{iS,cr,1,2}(1,:) ./ (CTs{iS,cr,1,2}(1,:) + CTs{iS,cr,1,1}(1,:)),'r');
plot(1:7, CTs{iS,cr,2,2}(2,:) ./ (CTs{iS,cr,2,2}(2,:) + CTs{iS,cr,2,1}(2,:)),'b:');
plot(1:7, CTs{iS,cr,2,2}(1,:) ./ (CTs{iS,cr,2,2}(1,:) + CTs{iS,cr,2,1}(1,:)),'r:');
hold on; hold all;
plot(1:7,propcw(1,:), 'r-h')
plot(1:7,propccw(1,:),'b-h')
plot(1:7,propcw(2,:), 'rh:')
plot(1:7,propccw(2,:),'bh:')


subplot(2,3,5);
hold  on; hold all;
plot(1:7,propcw(2,:), 'rh:')
plot(1:7,propccw(2,:),'bh:')
plot(1:7,PF{iS,cr,2}(2,:),'b-')
plot(1:7,PF{iS,cr,2}(1,:),'r-')
subplot(2,3,6);
hold on; hold all;
plot(1:7,propcw(1,:), 'r-h')
plot(1:7,propccw(1,:),'b-h')
plot(1:7,PF{iS,cr,1}(2,:),'b-')
plot(1:7,PF{iS,cr,1}(1,:),'r-')


