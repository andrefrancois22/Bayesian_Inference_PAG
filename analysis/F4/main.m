clear all; close all; clc;

% ==> directories
dataPath     = '/home/thomas/Desktop/UTAustin/Goris/pfc_code/pfc_data';
drc = '../../data/';
% ==> directory containing PF curve fit functions
pfc_functions_dr = '/home/thomas/Desktop/Bayesian_Inference_PAG/simulations/F1/pfc_functions/';
addpath(pfc_functions_dr)

% ==> dynamic range analysis types
dyn_types = {'popu','stim'};
% ==> estimate dynamic range per stimulus orientation? or per population?
dyn_type = dyn_types{1};

% ==> read in important metrics for final analysis
dvCatPerf = load([drc,'dvCatPerf.mat']);
dvCatPerf = dvCatPerf.dvCatPerf';
% ==> AIC Delta (from DV peak glm analysis in F2)
aic_d = load([drc,'aic_d.mat']);
aic_d = aic_d.aic_d;
% ==> log-likelihood ratios (from DV peak glm analysis in F2)
llr = load([drc,'llr.mat']);
llr = llr.llr';

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ==> PF curve estimation Set options
fitFlag         = 1;       % 0 or 1, where 1 means do the fit right now
nBootStraps     = 100;     % determines number of non parametric bootstraps if fitFlag == 1
options         = optimoptions('fmincon');
options.Display = 'off';
% Set bounds on model parameters
% Model fit to choice data split by task context (1 & 2: guess rate; 3: perceptual uncertainty; 4 & 5: decision criterion)
startVec_M1 = [0.04 0.04 1 0 0];
% startVec_V2 = [0.04 0.04 3 -15 15];
LB_M1(1,1)  = 0;                         UB_M1(1,1) = 0.05;      % lapse rate
LB_M1(2,1)  = 0;                         UB_M1(2,1) = 0.05;      % lapse rate
LB_M1(3,1)  = 0.1;                       UB_M1(3,1) = 10;        % perceptual uncertainty
LB_M1(4,1)  = -5;                        UB_M1(4,1) = 5;         % decision criterion (-10->10 range originally)
LB_M1(5,1)  = -5;                        UB_M1(5,1) = 5;         % decision criterion
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ==> trial counts and proportions in plots
tn  = nan(29,2,2,2,7);
tnp = nan(29,2,2,2,7);

% ==> delta bias and delta uncertainty metrics
db = nan(29,2);
dp = nan(29,2);

% ==> Psychometric functions (four per session)
% (<Session ID>_<contrast>_<dynamic range split>)
PF = cell([29,2,2]);
% ==> PF parameters
% (<Session ID>_<contrast>_<dynamic range split>)
PFp = cell([29,2,2]);

% ==> what is the session
for iS = 1:29 

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % dvs are model predicted DV trajectories [trial x time] (Category DVs)
    dvs = load([drc,'trial_DV_traj_iS_',num2str(iS),'.mat']);
    dvs = dvs.dv_cat;
    fprintf('Finished loading model predicted DV trajectories (Categorical DVs) for iS = %d... \n',iS)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ==> load session data
    if iS > 9
        load([dataPath,'/dataSet_',num2str(iS),'.mat']);
    elseif iS <= 9
        load([dataPath,'/dataSet_0',num2str(iS),'.mat']);
    end

    % ==> choice
    cho = S.beh.choiceCat;
    % ==> context, contrast, orientation indicator variables
    ctx = S.exp.taskContext; ctr = S.exp.stimContrast; ori = S.exp.stimOriDeg;
    % x axis is orientation (the values differ for FN and JP!). Use unique values
    or = sort(unique(ori),'ascend')'; cx = sort(unique(ctx),'ascend')'; cr = sort(unique(ctr),'ascend')'; bhs = sort(unique(cho),'ascend')';

    % ==> counts
    % (<choice>_<context>_<contrast>_<dynamic range split>_<orientation>)
    cnts = cell([2,2,2,2,2]);
    % ==> cell array containing indexing vectors
    % (<choice>_<context>_<contrast>_<dynamic range split>)
    idxs = cell([2,2,2,2]);

    % ==> option: split by dynamic range over entire population DVs
    if strcmp(dyn_type,'popu')
        % ==> dynamic range for session
        dynr = dynf(dvs); ldyn = (dynr(:,1) <  median(dynr)); hdyn = (dynr(:,1) >= median(dynr));
    end

    for bh = 1:length(bhs)             % --> behavioral choice (-1 or 1)
        for xc = 1:length(cx)          % --> prior context     (-1 or 1)
            for rc = 1:length(cr)      % --> stimulus contrast (min or max) (sorted in ascending order)
                for sp = 1:2           % --> dynamic range split (lo or hi) ldyn/hdyn                           

                    % ==> indexing vector for counts
                    idxs{bh}{xc}{rc}{sp} = ((ctx == cx(xc)) & (ctr == cr(rc)) & (cho == bhs(bh)));

                    % ==> stimulus orientation (7 values)
                    for th = 1:length(or)

                        % ==> option: split by dynamic range conditioned on identical stimulus values
                        if strcmp(dyn_type,'stim') 
                            % ==> dynamic range for session % => indexing variables - dynamic range
                            dv = dvs((ctx == cx(xc)) & (ctr == cr(rc)) & ori==or(th),:);
                            dynr = dynf(dv); ldyn=(dynf(dvs) <  median(dynr)); hdyn=(dynf(dvs) >= median(dynr));                           
                        end    
                        
                        % ==> dynamic range split indices (first column indexes low, second, high)
                        dri = [ldyn,hdyn];

                        % ==> tally counts
                        cnts{bh}{xc}{rc}{sp}{th} = sum(idxs{bh}{xc}{rc}{sp} & ori==or(th) & dri(:,sp));
                    end                    
                end
            end
        end    
    end
    % ==> fit psychometric functions for four PFs (by contrast, and DV dynamic range)
    for rc = 1:length(cr)
        for sp = 1:2      
            % ==> PF curve fit ~ input counts, orientations, initial parameter settings, options etc.
            [PF{iS}{rc}{sp}, PFp{iS}{rc}{sp}] = PF_fit_fun([cnts{1}{2}{rc}{sp}{:}], [cnts{1}{1}{rc}{sp}{:}], [cnts{2}{2}{rc}{sp}{:}], [cnts{2}{1}{rc}{sp}{:}], ...
             or, startVec_M1, LB_M1, UB_M1, options);                       
        end
        % ==> delta bias (use PF fit parameter vectors)
        db(iS,rc) = (PFp{iS}{rc}{1}(5) - PFp{iS}{rc}{1}(4)) - (PFp{iS}{rc}{2}(5) - PFp{iS}{rc}{2}(4));   
        % ==> delta perceptual uncertainty (use PF fit parameter vectors)
        dp(iS,rc)   = PFp{iS}{rc}{1}(3) - PFp{iS}{rc}{2}(3);         
    end
    % ==> update in the console
    fprintf('Computed psychometric functions for session %d of %d...\n',iS,29)                       
end

%%
close all;
% ==> Average delta bias by animal and by contrast 
% ==> correlation plots (delta bias and delta AIC) 
[~,~,~,~,~,~,~] = F4EFG(db, aic_d,'bias');

% ==> average delta perceptual uncertainty
% ==> correlation plots (delta perceptual uncertainty and delta AIC
[~,~,~,~,~,~,~] = F4EFG(dp, aic_d,'uncertainty');

%%
close all;
% ==> correlations between db and dp
F5(db, dp);

%%
% ==> example plot
% => Monkey colors
mclrs = {[1,0.75,0],[0.15,0.75,0.5]};
rgs = {1:13,14:29};
figure(); set(gcf,'color','white');
scatter(db(rgs{2},1), dp(rgs{2},1),75,'o','filled','markerfacecolor',mclrs{2},'markeredgecolor','w');
xlabel('Difference in decision bias');
ylabel('Difference in slope');
axis square;

%%
% ==> Useful supplementary plots
% => draw all PF curves and trial proportions for each session (4 PFs for each session)
SI_plots();

% % ==> ANCOVA
% dbz  = [db_lcr_z, db_hcr_z];
% dpz  = [dp_lcr_z, dp_hcr_z];
% aicz = [aic_d_z, aic_d_z];
% % => histogram to get grouping variable.
% [~,~,aicgv] = histcounts(aicz,8);
% % ==> run ANCOVA
% aoctool(dbz,dpz,aicgv)
