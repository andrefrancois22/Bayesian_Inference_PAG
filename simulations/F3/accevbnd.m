function [predPF_ddm, db, dp, prop_cw, prop_ccw, dvs_c_cw, dvs_i_cw, dvs_c_ccw, dvs_i_ccw, csensb_cw_ds, csensb_ccw_ds] = accevbnd(N, tm, bnd, sd, mus, fcs, ofs, M_FLAG, P_FLAG, N_FLAG, m, v, prf)

% ~~~~~~~~~~~~~~~ DDM EVIDENCE ACCUMULATION MODEL ~~~~~~~~~~~~~~~
% Runs the variants of the DDM (evidence accumulation model). Also produces
% simulated choices, median split analysis, and psychometric function
% fitting to the choices.
% ~~~~~~~~~~~~ INPUT ~~~~~~~~~~~~
% N              ==> Number of simulated trials
% tm             ==> Number of time points
% bnd            ==> Magnitude of the bound value. The bound value is only used when the model flag M_FLAG is equal to "BOUND"
% sd             ==> standard deviation for randn, which simulates the noise component of the momentary sensory evidence (the diffusion)
% mus            ==> simulated offsets for each of the 7 stimulus orientations
% fcs            ==> Prior offset factors (will influence drift rate due to integration)
% ofs            ==> The initial offset of the simulated DVs
% M_FLAG         ==> Evidence accumation type: BOUND or NO_BOUND (terminate at a bound, or terminate with the sign after tm timesteps)
% P_FLAG         ==> Model type: IMPULSE_PRIOR (only an initial offset), NO_IMPULSE (0 initial offset, followed by non-zero prior drift component), REG_DRIFT_PRIOR (an initial offset as well as a non-zero drift component following time 0)
% N_FLAG         ==> Trial noise parameter - cross-trial noise in the prior expectation
% m              ==> mean for log-normal cross-trial noise.
% v              ==> variance for log-normal cross-trial noise
% ~~~~~~~~~~~~ OUTPUT ~~~~~~~~~~~~
% predPF_ddm    ==> psychometric functions fit to the simulated choices following the median split analysis
% db            ==> Delta Bias results for each session: the difference in the decision criteria separating the low dynamic range DV PFs minus the decision criteria separating the high dynamic range DV PFs.
% dp            ==> Delta Uncertainty for each session: the difference in perceptual uncertainty between the PFs fit the low and the high dynamic range simulated DVs
% prop_cw       ==> choice proportions for simulated cw data
% prop_ccw      ==> choice proportions for the simulated ccw data
% dvs_c_cw      ==> simulated DVs for congruent cw cases
% dvs_c_ccw     ==> simulated DVs for congruent ccw cases
% dvs_i_cw      ==> simulated DVs for incongruent cw cases
% dvs_i_ccw     ==> simulated DVs for incongruent ccw cases
% csensb_cw_ds  ==> simulated DVs for all cw trials
% csensb_ccw_ds ==> simulated DVs for all ccw trials

% ==> delta bias
db = nan(1,length(fcs));
% ==> delta perceptual uncertainty
dp = nan(1,length(fcs));

% ==> store proportions for each dyn range split
prop_cw  = nan(length(fcs),3,7);
prop_ccw = nan(length(fcs),3,7);

% ==> counts (second index is context)
% <fcs><prior context><congruency><orientation><dynr condition>
cnts  = nan(length(fcs),2,2,7,2);

% ==> dynamic range split PF curve fits
%<fcs><dynamic range split case>
predPF_ddm = cell([length(fcs),2]);
% ==> PF parameters
paramEst_M1_ddm = cell([length(fcs),2]);

% ==> simulated DVs
csensb_cw_ds  = cell([length(fcs),7]);
csensb_ccw_ds = cell([length(fcs),7]);

% ==> dynamic range function
dynf = @(x) max([max(x, [], 2) - x(:,1), -(min(x, [], 2) - x(:,1))], [], 2);

% ==> orientation vector
% or_idx_v = [];

% ==> index
fci = 1;

for fc = fcs
        
    % ==> stimulus orientation
    for or = 1:7

        % ==> factor for prior offset in cw and ccw directions
        % ==> prior offsets
        pr_cw  =  ofs * fc; 
        pr_ccw = -ofs * fc; 

        % ==> parameter for changing prior strength
        pr_cw = pr_cw * prf;
        pr_ccw = pr_ccw * prf;
        
        if strcmp(P_FLAG,'REG_DRIFT_PRIOR')
            pr_cw_v  = repmat(pr_cw,  [tm,1]);  pr_cw_v(1)  =  ofs;
            pr_ccw_v = repmat(pr_ccw, [tm,1]);  pr_ccw_v(1) = -ofs;
            % % ==> set prior expectation step function for each trial
            pr_cw_vn  = repmat(pr_cw_v', [N,1]);  
            pr_ccw_vn = repmat(pr_ccw_v', [N,1]); 

            % ==> trial-by-trial noise?
            if strcmp(N_FLAG,'TRIAL_NOISE')
                mu = log((m^2)/sqrt(v+m^2));
                sigma = sqrt(log(v/(m^2)+1));

                epsi = lognrnd(mu,sigma,N,1);
                % ==> negative sign
                % epsi = epsi.*(-1).^randi(2,N,1);                
    
                % keyboard
                % ==> add single-trial noise to prior component
                pr_cw_vn  = pr_cw_vn  .* repmat(epsi,[1,tm]);
                pr_ccw_vn = pr_ccw_vn .* repmat(epsi,[1,tm]);   
                % keyboard
            end            
        end
        if strcmp(P_FLAG,'NO_IMPULSE')
            pr_cw_v  = repmat(pr_cw,  [tm,1]);  pr_cw_v(1)  = 0;
            pr_ccw_v = repmat(pr_ccw, [tm,1]);  pr_ccw_v(1) = 0;
            % % ==> set prior expectation step function for each trial
            pr_cw_vn  = repmat(pr_cw_v', [N,1]);  
            pr_ccw_vn = repmat(pr_ccw_v', [N,1]); 

            % ==> trial-by-trial noise?
            if strcmp(N_FLAG,'TRIAL_NOISE')
                mu = log((m^2)/sqrt(v+m^2));
                sigma = sqrt(log(v/(m^2)+1));

                epsi = lognrnd(mu,sigma,N,1);
                % ==> negative sign
                % epsi = epsi.*(-1).^randi(2,N,1);                
    
                % keyboard
                % ==> add single-trial noise to prior component
                pr_cw_vn  = pr_cw_vn  .* repmat(epsi,[1,tm]);
                pr_ccw_vn = pr_ccw_vn .* repmat(epsi,[1,tm]);   
                % keyboard
            end            
        end        
        % ==> case 'impulse function'
        if strcmp(P_FLAG,'IMPULSE_PRIOR')
            pr_cw_v  = zeros(tm,1);  pr_cw_v(1)  =  ofs * fc;
            pr_ccw_v = zeros(tm,1);  pr_ccw_v(1) = -ofs * fc; 
            % ==> set prior expectation step function for each trial
            pr_cw_vn  = repmat(pr_cw_v', [N,1]);
            pr_ccw_vn = repmat(pr_ccw_v', [N,1]);  

            % ==> trial-by-trial noise?
            if strcmp(N_FLAG,'TRIAL_NOISE')
                mu = log((m^2)/sqrt(v+m^2));
                sigma = sqrt(log(v/(m^2)+1));

                epsi = lognrnd(mu,sigma,N,1);
                % ==> negative sign
                % epsi = epsi.*(-1).^randi(2,N,1);                
    
                % keyboard
                % ==> add single-trial noise to prior component
                pr_cw_vn  = pr_cw_vn  .* repmat(epsi,[1,tm]);
                pr_ccw_vn = pr_ccw_vn .* repmat(epsi,[1,tm]);   
                % keyboard         
            end            
        end     

        % ==> randn
        rdn = randn(N,tm);
        % ==> Gaussian random walk
        sens_cw  = mus(or) + sd*rdn; %(momentary sens evidence)
        sens_ccw = mus(or) + sd*rdn;        
        
        % ==> cummulative sum over time (momentary evidence - )
        csensb_cw_d  = cumsum(sens_cw  + pr_cw_vn, 2);
        csensb_ccw_d = cumsum(sens_ccw + pr_ccw_vn,2);          

        % ==> once a bnd has been reached at time t for a given trial, set the remaining DV values
        % starting from t+1 to the bound for that trial
        % ==> cw cases
        csensb_cw_d = bndf(csensb_cw_d,bnd,N, M_FLAG);
        % ==> ccw context cases
        csensb_ccw_d = bndf(csensb_ccw_d,bnd,N, M_FLAG);

        if strcmp(M_FLAG,'BOUND')
            %==> simulating responses and dynamic range analysis
            % ==> responses for bound crossings
            % ==> cw
            cw_dv_cw_r =    (csensb_cw_d(:,end)  ==   bnd);
            cw_dv_ccw_r =  -(csensb_cw_d(:,end)  ==  -bnd);
            % ==> ccw
            ccw_dv_cw_r =   (csensb_ccw_d(:,end) ==  bnd);
            ccw_dv_ccw_r = -(csensb_ccw_d(:,end) == -bnd);
            
        elseif strcmp(M_FLAG,'NO_BOUND')   
            % ==> no bound - responses are sign at the end of DV
            % ==> cw
            cw_dv_cw_r =    (sign(csensb_cw_d(:,end))  ==   1);
            cw_dv_ccw_r =  -(sign(csensb_cw_d(:,end))  ==  -1);
            % ==> ccw
            ccw_dv_cw_r =   (sign(csensb_ccw_d(:,end)) ==  1);
            ccw_dv_ccw_r = -(sign(csensb_ccw_d(:,end)) == -1);                        
        end   
        
        % ==> some sanity checks
        assert(sum([cw_dv_cw_r  & cw_dv_ccw_r]) == 0);
        assert(sum([ccw_dv_cw_r & ccw_dv_ccw_r]) == 0);

        % ==> create 2AFC response vector (for cw context DVs)
        cw_r  = [cw_dv_cw_r  + cw_dv_ccw_r]; 
        % ==> for ccw context DVs
        ccw_r = [ccw_dv_cw_r + ccw_dv_ccw_r];

        % ==> handle DVs that never reached a bound...
        % ==> use sign of DV
        cw_r(cw_r == 0)   = sign(csensb_cw_d(cw_r==0,end));
        ccw_r(ccw_r == 0) = sign(csensb_ccw_d(ccw_r==0,end));            

        % ==> incongruent trials (cw context yields a ccw response, e.g. '-1' )
        dvs_i_cw = csensb_cw_d(cw_r==-1,:);
        % ==> congruent trials (cw context)
        dvs_c_cw = csensb_cw_d(cw_r== 1,:);
 
        % ==> incongruent trials (ccw context yields a cw response, e.g. '-1' )
        dvs_i_ccw = csensb_ccw_d(ccw_r== 1,:);
        % ==> congruent trials (ccw context yields a ccw response e.g. '-1')
        dvs_c_ccw = csensb_ccw_d(ccw_r==-1,:);

        % ==> cw context, and stim orientation dynr median
        med_cw  = median(dynf(csensb_cw_d));
        % ==> ccw context, and stim orientation dynr median
        med_ccw = median(dynf(csensb_ccw_d));

        % ==> proportions => high dynamic range DVs
        prop_cw(fci,2,or)  = sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) > med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) > med_cw) == -1));
        prop_ccw(fci,2,or) = sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) > med_ccw) == -1));
        % ==> proportions => low dynamic range DVs
        prop_cw(fci,1,or)  = sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)     / (sum(cw_r(dynf(csensb_cw_d) <= med_cw) == 1)    + sum(cw_r(dynf(csensb_cw_d) <= med_cw) == -1));
        prop_ccw(fci,1,or) = sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1)  / (sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == 1) + sum(ccw_r(dynf(csensb_ccw_d) <= med_ccw) == -1));
        % ==> proportions without split
        prop_cw(fci,3,or)  = sum(cw_r  == 1) / (sum(cw_r  == 1) + sum(cw_r  == -1));
        prop_ccw(fci,3,or) = sum(ccw_r == 1) / (sum(ccw_r == 1) + sum(ccw_r == -1));

        % ==> counts ( ==> high dynr)
        % ==> congruent cases
        cnts(fci,1,1,or,1) = sum(cw_dv_cw_r(  dynf(csensb_cw_d)  > med_cw)  ==  1);
        cnts(fci,2,1,or,1) = sum(ccw_dv_ccw_r(dynf(csensb_ccw_d) > med_ccw) == -1);
        % ==> incongruent cases
        cnts(fci,1,2,or,1) = sum(cw_dv_ccw_r(dynf(csensb_cw_d)  > med_cw)  == -1);
        cnts(fci,2,2,or,1) = sum(ccw_dv_cw_r(dynf(csensb_ccw_d) > med_ccw) ==  1);        
        % ==> counts ( ==> low dynr)
        % ==> congruent cases
        cnts(fci,1,1,or,2) = sum(cw_dv_cw_r(  dynf(csensb_cw_d)  <= med_cw)  ==  1);
        cnts(fci,2,1,or,2) = sum(ccw_dv_ccw_r(dynf(csensb_ccw_d) <= med_ccw) == -1);
        % ==> incongruent cases
        cnts(fci,1,2,or,2) = sum(cw_dv_ccw_r(dynf(csensb_cw_d)  <= med_cw)  == -1);
        cnts(fci,2,2,or,2) = sum(ccw_dv_cw_r(dynf(csensb_ccw_d) <= med_ccw) ==  1);        

        % ==> save all simulated DVs
        csensb_cw_ds{fci,or}  = csensb_cw_d;
        csensb_ccw_ds{fci,or} = csensb_ccw_d;

        % ==> will use this to compute F3A scatterplots filtering out all
        % oriented cases (keeping only simulated neutral orientations)
        %or_idx_v = [or_idx_v; repmat(or,[N,1])];
    end
    
    % ==> for dynr case (1 - high or 2 - low)
    for dri = 1:2

        % ==> fit psychometric functions
        % ==> use identical db and dp metrics (must do this for dynamic range proportions too)
        fprintf('computing psychometric functions...\n')
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        options         = optimoptions('fmincon');
        options.Display = 'off';
        % Set bounds on model parameters
        % Model fit to choice data split by task context (1 & 2: guess rate; 3: perceptual uncertainty; 4 & 5: decision criterion)
        startVec_M1 = [0.04 0.04 1 0 0]; 
    
        LB_M1(1,1)  = 0;                         UB_M1(1,1) = 0.05;      % lapse rate
        LB_M1(2,1)  = 0;                         UB_M1(2,1) = 0.05;      % lapse rate
        LB_M1(3,1)  = 0.0;                       UB_M1(3,1) = 5;        % perceptual uncertainty
        LB_M1(4,1)  = -10;                       UB_M1(4,1) = 10;         % decision criterion (-10->10 range originally)
        LB_M1(5,1)  = -10;                       UB_M1(5,1) = 10;         % decision criterion
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % ==> context 1 ccw response (incongruent) & ==> context 2 ccw response (congruent)
        nccw = [squeeze(cnts(fci,1,2,:,dri))'; squeeze(cnts(fci,2,1,:,dri))']; 
        % ==> context 1 cw response (congruent)    & ==> context 2 cw response (incongruent)
        ncw  = [squeeze(cnts(fci,1,1,:,dri))'; squeeze(cnts(fci,2,2,:,dri))']; 
        % ==> pass through PF function fit
        obFun    = @(paramVec) giveNLL(paramVec, mus, nccw, ncw, 'M1');   

        NLL = inf;
        while NLL == inf
            % ==> run optimization
            paramEst_M1_ddm{fci,dri} = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);        
            % ==> return pf
            [NLL, predPF_ddm{fci,dri}] = giveNLL(paramEst_M1_ddm{fci,dri}, mus, nccw, ncw, 'M1');
            if NLL > 250 % ==> somewhat arbitrary limit, but saves us from local minima
                fprintf(['NLL = ',num2str(NLL),', run fmincon again for better results...\n'])
                % ==> NLL bigger than 100 usually means fmincon converged
                % to a local optimum and PF curve fitting will not be good
                % ==> don't break out of while loop - try again
                NLL = inf;                
            end
            % ==> try optimization again with slightly jittered initial
            % parameter settings by adding a small epsilon of Gaussian
            % noise
            startVec_M1 = startVec_M1 + randn(1,5);
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
    end

    % ==> store bias estimates (map) - calculating as:
    % the distance between the points of subjective equality 
    % for the cw and ccw prior context PF curves
    db(fci) = (paramEst_M1_ddm{fci,2}(end) - paramEst_M1_ddm{fci,2}(end-1)) - (paramEst_M1_ddm{fci,1}(end) - paramEst_M1_ddm{fci,1}(end-1));   
    % ==> perceptual uncertainty
    dp(fci) = paramEst_M1_ddm{fci,2}(3) - paramEst_M1_ddm{fci,1}(3);   

    % keyboard

    % ==> update
    fprintf('Simulated PFs for drift rate and offset parameter %d of %d...\n', fci, length(fcs));
    % ==> increment
    fci = fci + 1;    
end
