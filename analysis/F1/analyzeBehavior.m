% Start with clean slate
clearvars
clc

% Set Paths
thisPath     = pwd; %strcat('/Users/g4344/Desktop/Research/Finished/Charlton & Goris (2024)/pfc_analysis');
dataPath     = strcat('../../data/pfc_data/');
functionPath = strcat('../../simulations/F1/pfc_functions/');

% Choose example data-set
plotExpPrior = [7, 24];  % F: 1-13; JP: 14-29
plotExpRule  = [8, 24];  % F: 1-13; JP: 14-29

% Set options
fitFlag         = 0;       % 0 or 1, where 1 means do the fit right now
nBootStraps     = 100;     % determines number of non parametric bootstraps if fitFlag == 1
options         = optimoptions('fmincon');
options.Display = 'off';

% Set bounds on model parameters
% Model fit to choice data split by task context (1 & 2: guess rate; 3: perceptual uncertainty; 4 & 5: decision criterion)
startVec_M1 = [0.01 0.01 1 0 0];
LB_M1(1,1)  = 0;                        UB_M1(1,1) = 0.05;      % lapse rate
LB_M1(2,1)  = 0;                        UB_M1(2,1) = 0.05;      % lapse rate
LB_M1(3,1)  = 0.1;                      UB_M1(3,1) = 10;        % perceptual uncertainty
LB_M1(4,1)  = -10;                      UB_M1(4,1) = 10;        % decision criterion
LB_M1(5,1)  = -10;                      UB_M1(5,1) = 10;        % decision criterion

% Model fit to choice data split by mapping rule (1 & 2: guess rate; 3 & 4: perceptual uncertainty; 5 & 6: decision criterion)
startVec_M2 = [0.01 0.01 1 1 0 0];
LB_M2(1,1)  = 0;                        UB_M2(1,1) = 0.05;      % lapse rate
LB_M2(2,1)  = 0;                        UB_M2(2,1) = 0.05;      % lapse rate
LB_M2(3,1)  = 0.1;                      UB_M2(3,1) = 10;        % perceptual uncertainty
LB_M2(4,1)  = 0.1;                      UB_M2(4,1) = 10;        % perceptual uncertainty
LB_M2(5,1)  = -10;                      UB_M2(5,1) = 10;        % decision criterion
LB_M2(6,1)  = -10;                      UB_M2(6,1) = 10;        % decision criterion


% Make list of file names
cd(dataPath)
fileNameList = dir('dataSet_*.mat');


% Loop through recording sessions
for iS = 1:numel(fileNameList)
    
    % Specify load name
    loadName = fileNameList(iS).name;
    
    % Print some output
    fprintf('Analyzing data set %d of %d... \n', iS, numel(fileNameList))
        
    % Load raw data
    cd(dataPath)
    load(loadName);
    
    if fitFlag == 1
        
        % Identify variables to select specific choice data (structure: 2 priors X 2 rules X 2 contrasts x 7 orientations)
        contrastValues = unique(S.exp.stimContrast);
        priorValues    = unique(S.exp.taskContext);
        oriValues      = unique(S.exp.stimOriDeg);
        ruleValues     = unique(S.exp.respContext);
        
        % Compute the proportion of saccades that expresses a "CW" choice
        propRightCW = sum(S.exp.respContext == ruleValues(1) & S.beh.choiceCat == 1)/sum(S.exp.respContext == ruleValues(1));
        propLeftCW  = sum(S.exp.respContext == ruleValues(2) & S.beh.choiceCat == 1)/sum(S.exp.respContext == ruleValues(2));

        % Compute the average response time under each mapping rule
        respTimeRule1 = nanmean(S.beh.respTime(S.exp.respContext == ruleValues(1)) - S.exp.Report_on(S.exp.respContext == ruleValues(1)));
        respTimeRule2 = nanmean(S.beh.respTime(S.exp.respContext == ruleValues(2)) - S.exp.Report_on(S.exp.respContext == ruleValues(2)));

        % Compute performance on the "catch" trials under each mapping rule
        perfCatchRule1 = mean(S.beh.choiceEval(S.exp.respContext == ruleValues(1) & S.exp.stimContrast == contrastValues(2) & abs(S.exp.stimOriDeg) == oriValues(end)));
        perfCatchRule2 = mean(S.beh.choiceEval(S.exp.respContext == ruleValues(2) & S.exp.stimContrast == contrastValues(2) & abs(S.exp.stimOriDeg) == oriValues(end)));

        % Compute performance on all trials under each mapping rule
        perfAllRule1 = mean(S.beh.choiceEval(S.exp.respContext == ruleValues(1)));
        perfAllRule2 = mean(S.beh.choiceEval(S.exp.respContext == ruleValues(2)));

        for iC = 1:numel(contrastValues)
                        
            for iP = 1:numel(priorValues) % also used for mapping rules

                % Create useful trial indices
                indRule  = (S.exp.respContext == ruleValues(iP));
                indPrior = (S.exp.taskContext == priorValues(iP));
                
                for iO = 1:numel(oriValues)
                                        
                    % Create useful trial index
                    indStim = (S.exp.stimContrast == contrastValues(iC) & S.exp.stimOriDeg == oriValues(iO));                                       
                    
                    % Split out choice data by response context (mapping rule 1 vs 2)
                    nCCW_rule(iP,iO) = sum(S.beh.choiceCat(indStim & indRule) == -1);
                    nCW_rule(iP,iO)  = sum(S.beh.choiceCat(indStim & indRule) == 1);
                  
                    % Split out choice data by task context (CW vs CCW condition)
                    nCCW_prior(iP,iO) = sum(S.beh.choiceCat(indStim & indPrior) == -1);
                    nCW_prior(iP,iO)  = sum(S.beh.choiceCat(indStim & indPrior) == 1);
                    nLW_prior(iP,iO)  = sum(S.beh.choiceDir(indStim & indPrior) == -1);
                    nRW_prior(iP,iO)  = sum(S.beh.choiceDir(indStim & indPrior) == 1);
                    
                    % Characterize task design 
                    nStimCCW_prior(iP,iO) = sum(S.exp.stimOriDeg(indStim & indPrior) < 0) + .5*sum(S.exp.stimOriDeg(indStim & indPrior) == 0);
                    nStimCW_prior(iP,iO)  = sum(S.exp.stimOriDeg(indStim & indPrior) > 0) + .5*sum(S.exp.stimOriDeg(indStim & indPrior) == 0);
                    nStimLW_prior(iP,iO)  = sum(((S.exp.respContext(indStim & indPrior) == 0) & (S.exp.stimOriDeg(indStim & indPrior) < 0))|((S.exp.respContext(indStim & indPrior) ~= 0) & (S.exp.stimOriDeg(indStim & indPrior) > 0))) + .5*sum(S.exp.stimOriDeg(indStim & indPrior) == 0);
                    nStimRW_prior(iP,iO)  = sum(((S.exp.respContext(indStim & indPrior) == 1) & (S.exp.stimOriDeg(indStim & indPrior) < 0))|((S.exp.respContext(indStim & indPrior) ~= 1) & (S.exp.stimOriDeg(indStim & indPrior) > 0))) + .5*sum(S.exp.stimOriDeg(indStim & indPrior) == 0);                                       
                end
            end
                       
            % Summarize ideal response distributions
            taskDirLogOdds_prior{iC} = log(sum(nStimLW_prior, 2)./sum(nStimRW_prior, 2));
            taskCatLogOdds_prior{iC} = log(sum(nStimCCW_prior, 2)./sum(nStimCW_prior, 2));
            
            % Summarize choice behavior for ambiguous stimuli
            dirLogOdds_prior{iC} = log(nLW_prior(:, 4)./nRW_prior(:, 4));
            catLogOdds_prior{iC} = log(nCCW_prior(:, 4)./nCW_prior(:, 4));

            % Choice data split by task context
            propCW_prior{iC}  = nCW_prior./(nCW_prior + nCCW_prior);
            nTrials_prior{iC} = nCW_prior + nCCW_prior;

            % Choice data split by response rule
            propCW_rule{iC}  = nCW_rule./(nCW_rule + nCCW_rule);
            nTrials_rule{iC} = nCW_rule + nCCW_rule;
                        
            % Fit model to estimate sensitivity and bias
            cd(functionPath)
            
            obFun               = @(paramVec) giveNLL(paramVec, oriValues, nCCW_prior, nCW_prior, 'M1');         
            paramEst_M1{iS}{iC} = fmincon(obFun, startVec_M1, [], [], [], [], LB_M1, UB_M1, [], options);
            
            obFun               = @(paramVec) giveNLL(paramVec, oriValues, nCCW_rule, nCW_rule, 'M2');
            paramEst_M2{iS}{iC} = fmincon(obFun, startVec_M2, [], [], [], [], LB_M2, UB_M2, [], options);
                        
            % Fit boostrapped data to get confidence intervals
            for iBS = 1:nBootStraps
                
                % Bootstrap by drawing from binomial distribution
                nCW_prior_bs  = binornd(nCCW_prior + nCW_prior, nCW_prior./(nCCW_prior + nCW_prior));
                nCCW_prior_bs = (nCCW_prior + nCW_prior) - nCW_prior_bs;

                nCW_rule_bs  = binornd(nCCW_rule + nCW_rule, nCW_rule./(nCCW_rule + nCW_rule));
                nCCW_rule_bs = (nCCW_rule + nCW_rule) - nCW_rule_bs;
                               
                obFun           = @(paramVec) giveNLL(paramVec, oriValues, nCCW_prior_bs, nCW_prior_bs, 'M1');
                paramsBS_M1{iC} = fmincon(obFun, paramEst_M1{iS}{iC}, [], [], [], [], LB_M1, UB_M1, [], options);
                
                obFun           = @(paramVec) giveNLL(paramVec, oriValues, nCCW_rule_bs, nCW_rule_bs, 'M2');
                paramsBS_M2{iC} = fmincon(obFun, paramEst_M2{iS}{iC}, [], [], [], [], LB_M2, UB_M2, [], options);
                                
                % Document outcome analysis
                S.beh.model(1).oriSens.bs{iC}(iBS)   = 1./paramsBS_M1{iC}(3);
                S.beh.model(1).oriBias.bs{iC}(iBS)   = diff(paramsBS_M1{iC}([5,4]));                
                S.beh.model(2).oriSens.bs{iC}(iBS,:) = 1./paramsBS_M2{iC}([3,4]);
                S.beh.model(2).oriPSE.bs{iC}(iBS,:)  = paramsBS_M2{iC}([5,6]);
            end            
        end
                       
        % Document outcome analysis
        S.beh.respTimeRule1        = respTimeRule1;
        S.beh.respTimeRule2        = respTimeRule2;
        S.beh.propRightCW          = propRightCW;
        S.beh.propLeftCW           = propLeftCW;
        S.beh.perfCatchRule1       = perfCatchRule1;
        S.beh.perfCatchRule2       = perfCatchRule2;
        S.beh.perfAllRule1         = perfAllRule1;
        S.beh.perfAllRule2         = perfAllRule2;
        
        S.beh.nTrials_prior        = nTrials_prior;
        S.beh.nTrials_rule         = nTrials_rule;
        S.beh.propCW_prior         = propCW_prior;
        S.beh.propCW_rule          = propCW_rule;
        
        S.beh.dirLogOdds           = dirLogOdds_prior;
        S.beh.catLogOdds           = catLogOdds_prior;
        S.exp.taskDirLogOdds       = taskDirLogOdds_prior;
        S.exp.taskCatLogOdds       = taskCatLogOdds_prior;
        
        S.beh.model(1).comment     = 'Two sets of five parameters (one per contrast) describe the proportion of CW choices as a function of stimulus orientation and prior context';
        S.beh.model(1).params      = paramEst_M1{iS};
        S.beh.model(1).oriSens.est = [1./paramEst_M1{iS}{1}(3), 1./paramEst_M1{iS}{2}(3)];
        S.beh.model(1).oriBias.est = [diff(paramEst_M1{iS}{1}([5,4])), diff(paramEst_M1{iS}{2}([5,4]))];

        S.beh.model(2).comment     = 'Two sets of six parameters (one per contrast) describe the proportion of CW choices as a function of stimulus orientation and mapping rule';
        S.beh.model(2).params      = paramEst_M2{iS};
        S.beh.model(2).oriSens.est = [1./paramEst_M2{iS}{1}([3,4]); 1./paramEst_M2{iS}{2}([3,4])];
        S.beh.model(2).oriPSE.est  = [paramEst_M2{iS}{1}([5,6]); paramEst_M2{iS}{2}([5,6])];
        
        % Save analysis
        cd(dataPath)
        save(loadName, 'S');
    end
    
    % Computations for plotting
    dataF(iS)  = eval(S.general.expDate(1:2)) < 21;   % 2020 experiments (Friedrich)
    dataJP(iS) = eval(S.general.expDate(1:2)) >= 21;  % 2021 and 2022 experiments (JP)

    % The proportion of saccades that express "CW" choice
    propRightCW(iS) = S.beh.propRightCW;
    propLeftCW(iS)  = S.beh.propLeftCW;
    
    % The response time under each mapping rule
    respTimeRule1(iS) = S.beh.respTimeRule1;
    respTimeRule2(iS) = S.beh.respTimeRule2;
    
    % Performance on the catch trials
    perfCatchRule1(iS) = S.beh.perfCatchRule1;
    perfCatchRule2(iS) = S.beh.perfCatchRule2;

    % % Performance on all trials
    % perfAllRule1(iS) = S.beh.perfAllRule1;
    % perfAllRule2(iS) = S.beh.perfAllRule2;

    % Orientation sensitivity per rule
    oriSensLC_rule(iS,:)       = S.beh.model(2).oriSens.est(1,:);
    oriSensHC_rule(iS,:)       = S.beh.model(2).oriSens.est(2,:);
    oriSensLC_IQR_rule(iS,:,:) = prctile(S.beh.model(2).oriSens.bs{1}, [25 75]);
    oriSensHC_IQR_rule(iS,:,:) = prctile(S.beh.model(2).oriSens.bs{2}, [25 75]);

    % Orientation PSE per rule
    oriPSELC_rule(iS,:) = S.beh.model(2).oriPSE.est(1,:);
    oriPSEHC_rule(iS,:) = S.beh.model(2).oriPSE.est(2,:);
    
    % Ideal response distribution given task statistics
    taskDirLogOddsLC(iS,:) = S.exp.taskDirLogOdds{1};
    taskDirLogOddsHC(iS,:) = S.exp.taskDirLogOdds{2};
    taskCatLogOddsLC(iS,:) = S.exp.taskCatLogOdds{1};
    taskCatLogOddsHC(iS,:) = S.exp.taskCatLogOdds{2};
    
    % Observed response distribution
    dirLogOddsLC(iS,:) = S.beh.dirLogOdds{1};
    dirLogOddsHC(iS,:) = S.beh.dirLogOdds{2};
    catLogOddsLC(iS,:) = S.beh.catLogOdds{1};
    catLogOddsHC(iS,:) = S.beh.catLogOdds{2};
    
    % Number of completed trials
    trialsDoneCCWLC(iS,:) = S.beh.nTrials_prior{1}(1,:);
    trialsDoneCWLC(iS,:)  = S.beh.nTrials_prior{1}(2,:);
    trialsDoneCCWHC(iS,:) = S.beh.nTrials_prior{2}(1,:);
    trialsDoneCWHC(iS,:)  = S.beh.nTrials_prior{2}(2,:);
    
    % Model based estimate of orientation sensitivity and bias
    oriSensLC(iS) = S.beh.model(1).oriSens.est(1);
    oriSensHC(iS) = S.beh.model(1).oriSens.est(2);
    oriBiasLC(iS) = S.beh.model(1).oriBias.est(1);
    oriBiasHC(iS) = S.beh.model(1).oriBias.est(2);
    
    oriSensLC_IQR(iS,:) = prctile(S.beh.model(1).oriSens.bs{1}, [25 75]);
    oriSensHC_IQR(iS,:) = prctile(S.beh.model(1).oriSens.bs{2}, [25 75]);
    oriBiasLC_IQR(iS,:) = prctile(S.beh.model(1).oriBias.bs{1}, [25 75]);
    oriBiasHC_IQR(iS,:) = prctile(S.beh.model(1).oriBias.bs{2}, [25 75]);
end

%% Compute some summary statistics and conduct some statistical tests
catBiasDiff_F  = median([catLogOddsLC(dataF,1); catLogOddsHC(dataF,1)] - [catLogOddsLC(dataF,2); catLogOddsHC(dataF,2)]);
catBiasDiff_JP = median([catLogOddsLC(dataJP,1); catLogOddsHC(dataJP,1)] - [catLogOddsLC(dataJP,2); catLogOddsHC(dataJP,2)]);
dirBiasDiff_F  = median([dirLogOddsLC(dataF,1); dirLogOddsHC(dataF,1)] - [dirLogOddsLC(dataF,2); dirLogOddsHC(dataF,2)]);
dirBiasDiff_JP = median([dirLogOddsLC(dataJP,1); dirLogOddsHC(dataJP,1)] - [dirLogOddsLC(dataJP,2); dirLogOddsHC(dataJP,2)]);

pCatBiasDiff_F  = signrank([catLogOddsLC(dataF,1); catLogOddsHC(dataF,1)] - [catLogOddsLC(dataF,2); catLogOddsHC(dataF,2)]);
pCatBiasDiff_JP = signrank([catLogOddsLC(dataJP,1); catLogOddsHC(dataJP,1)] - [catLogOddsLC(dataJP,2); catLogOddsHC(dataJP,2)]);
pDirBiasDiff_F  = signrank([dirLogOddsLC(dataF,1); dirLogOddsHC(dataF,1)] - [dirLogOddsLC(dataF,2); dirLogOddsHC(dataF,2)]);
pDirBiasDiff_JP = signrank([dirLogOddsLC(dataJP,1); dirLogOddsHC(dataJP,1)] - [dirLogOddsLC(dataJP,2); dirLogOddsHC(dataJP,2)]);


%% Make figures
for iA = 1:2       % Two animals 
    
    if iA == 1
        dataSets   = dataF;
        monkeyName = 'Friedrich';
    elseif iA == 2
        dataSets   = dataJP;
        monkeyName = 'JP';
    end

    % Compute some monkey-specific statistics
    [rSensBiasAll, pSensBiasAll] = corr(log2([oriSensLC(dataSets), oriSensHC(dataSets)]'), [oriBiasLC(dataSets), oriBiasHC(dataSets)]', 'type', 'Spearman');
    [rSensBiasLC, pSensBiasLC]   = corr(log2(oriSensLC(dataSets)'), oriBiasLC(dataSets)', 'type', 'Spearman');
    [rSensBiasHC, pSensBiasHC]   = corr(log2(oriSensHC(dataSets)'), oriBiasHC(dataSets)', 'type', 'Spearman');
    
    [pOriSenDiff, ~] = signrank(1 - ([oriSensLC_rule(dataSets,2); oriSensHC_rule(dataSets,2)]./[oriSensLC_rule(dataSets,1); oriSensHC_rule(dataSets,1)]));
    
    % Set first figure
    set(figure(2*(iA-1)+1), 'OuterPosition', [100 100 2000 1000])
    
    subplot(2,4,1)
    plot([0 3], [.5 .5], 'k--', 'linewidth', 1)
    hold on, box off, axis square
    plot(1 + .1*randn(sum(dataSets), 1), propRightCW(dataSets), 'ko', 'markersize', 9, 'markerfacecolor', [.5 .5 .5])
    plot(2 + .1*randn(sum(dataSets), 1), propLeftCW(dataSets), 'ko', 'markersize', 9, 'markerfacecolor', [.5 .5 .5])
    axis([0 3 0 1])
    title(strcat("Choice behavior ", monkeyName))
    xlabel("Mapping rule")
    ylabel("Clockwise choices (%)")

    subplot(2,4,2)
    plot([0.1 0.25], [0.1 0.25], 'k--')
    hold on, box off, axis square
    plot(respTimeRule1(dataSets), respTimeRule2(dataSets), 'ko', 'markersize', 9, 'markerfacecolor', [.5 .5 .5])
    axis([0.1 0.25 0.1 0.25])
    xlabel("Response speed rule 1 (sec)")
    ylabel("Response speed rule 2 (sec)") 
    
    subplot(2,4,3)
    plot(perfCatchRule1(dataSets), perfCatchRule2(dataSets), 'ko', 'markersize', 9, 'markerfacecolor', [.5 .5 .5])
    hold on, box off, axis square
    plot([.5 .5], [.45 1], 'k-')
    plot([.45 1], [.5 .5], 'k-')
    plot([.45 1], [.45 1], 'k--')
    axis([.45 1 .45 1])
    xlabel("Catch trials rule 1 (% correct)")
    ylabel("Catch trials rule 2 (% correct)")
    
    subplot(2,4,6)
    plot(log2(oriSensLC_rule(dataSets,1)), log2(oriSensLC_rule(dataSets,2)), 'ko', 'markersize', 9, 'markerfacecolor', [.85 .85 .85])
    hold on, box off, axis square
    plot(log2(oriSensHC_rule(dataSets,1)), log2(oriSensHC_rule(dataSets,2)), 'ko', 'markersize', 9, 'markerfacecolor', [0 0 0])
    plot(log2(oriSensLC_IQR_rule(dataSets,:,1))', [log2(oriSensLC_rule(dataSets,2)), log2(oriSensLC_rule(dataSets,2))]', 'k-')
    plot([log2(oriSensLC_rule(dataSets,1)), log2(oriSensLC_rule(dataSets,1))]', log2(oriSensLC_IQR_rule(dataSets,:,2))', 'k-') 
    plot(log2(oriSensHC_IQR_rule(dataSets,:,1))', [log2(oriSensHC_rule(dataSets,2)), log2(oriSensHC_rule(dataSets,2))]', 'k-')
    plot([log2(oriSensHC_rule(dataSets,1)), log2(oriSensHC_rule(dataSets,1))]', log2(oriSensHC_IQR_rule(dataSets,:,2))', 'k-')
    plot([log2(.25) log2(4)], [log2(.25) log2(4)], 'k--')
    axis([log2(.25) log2(4) log2(.25) log2(4)])
    legend('Low contrast', 'High contrast', 'IQR', 'location', 'NorthWest')
    xlabel("Orientation sensitivity rule 1 (log_2 basis)")
    ylabel("Orientation sensitivity rule 2(log_2 basis)")
    title(sprintf('P = %0.2g (Wilcoxson)', pOriSenDiff))

    subplot(2,4,7)
    plot(oriPSELC_rule(dataSets,1), oriPSELC_rule(dataSets,2), 'ko', 'markersize', 9, 'markerfacecolor', [.85 .85 .85])
    hold on, box off, axis square
    plot(oriPSEHC_rule(dataSets,1), oriPSEHC_rule(dataSets,2), 'ko', 'markersize', 9, 'markerfacecolor', [0 0 0])
    plot([-3 3], [-3 3], 'k--')
    axis([-3 3 -3 3])
    legend('Low contrast', 'High contrast', 'location', 'NorthWest')
    xlabel("Orientation PSE rule 1 (log_2 basis)")
    ylabel("Orientation PSE rule 2(log_2 basis)")
       
    
    % Plot example psychometric functions
    cd(dataPath)
    load(fileNameList(plotExpRule(iA)).name);
    cd(functionPath)
        
    % Computations for plotting
    nSamples    = 100;
    stimPlot    = linspace(min(S.general.oriValuesDeg), max(S.general.oriValuesDeg), nSamples);
    [~, predPF] = giveNLL(S.beh.model(2).params{2}, stimPlot, zeros(2,nSamples), zeros(2,nSamples), 'M2');
    
    % Plot example dataset
    subplot(2,4,5)
    plot(stimPlot, predPF(1,:), '-', 'linewidth', 2, 'color', [1 .5 0])
    hold on, box off, axis square
    plot(stimPlot, predPF(2,:), '-', 'linewidth', 2, 'color', [0 .7 .2])
    plot([-5 5], [.5 .5], 'k--')
    plot([0 0], [0 1], 'k--')
    axis([-5 5 0 1])
    xlabel("Stimulus orientation (deg)")
    ylabel("Proportion clockwise")
    title("High stimulus contrast")
    
    for iO = 1:numel(S.general.oriValuesDeg)
        plot(S.general.oriValuesDeg(iO), S.beh.propCW_rule{2}(1,iO), 'ko', 'markerfacecolor', [1 .5 0], 'markersize', round(S.beh.nTrials_rule{2}(1,iO)/5)+1)
        plot(S.general.oriValuesDeg(iO), S.beh.propCW_rule{2}(2,iO), 'ko', 'markerfacecolor', [0 .7 .2], 'markersize', round(S.beh.nTrials_rule{2}(2,iO)/5)+1)
    end
    legend('Rule 1', 'Rule 2', 'location', 'NorthWest')
    cd(thisPath)

    
    % Set second figure 
    set(figure(2*(iA-1)+2), 'OuterPosition', [100 100 2000 1000])
    
    subplot(2,4,1)
    plot(S.general.oriValuesDeg, .5*mean(trialsDoneCCWLC(dataSets,:) + trialsDoneCCWHC(dataSets,:)), '-', 'linewidth', 3, 'color', [1 0 0])
    hold on, box off, axis square
    plot(S.general.oriValuesDeg, .5*mean(trialsDoneCWLC(dataSets,:)  + trialsDoneCWHC(dataSets,:)), '-', 'linewidth', 3, 'color', [0 .7 .9])
    plot(S.general.oriValuesDeg, trialsDoneCCWLC(dataSets,:), 'ko', 'markersize', 8, 'markerfacecolor', [1 .75 .75])
    plot(S.general.oriValuesDeg, trialsDoneCCWHC(dataSets,:), 'ko', 'markersize', 8, 'markerfacecolor', [1 0 0])
    plot(S.general.oriValuesDeg, trialsDoneCWLC(dataSets,:), 'ko', 'markersize', 8, 'markerfacecolor', [.5 .9 1])
    plot(S.general.oriValuesDeg, trialsDoneCWHC(dataSets,:), 'ko', 'markersize', 8, 'markerfacecolor', [0 .7 .9])
    axis([-5 5 0 350])
    legend('Context 1', 'Context 2', 'location', 'NorthWest')
    title(strcat("Choice behavior ", monkeyName))
    xlabel("Stimulus orientation (deg)")
    ylabel("Completed trials")
    
    subplot(2,4,2)
    plot(taskCatLogOddsLC(dataSets,1), taskDirLogOddsLC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [1 .75 .75])
    hold on, box off, axis square
    plot(taskCatLogOddsHC(dataSets,1), taskDirLogOddsHC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [1 0 0])
    plot(taskCatLogOddsLC(dataSets,2), taskDirLogOddsLC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [.5 .9 1])
    plot(taskCatLogOddsHC(dataSets,2), taskDirLogOddsHC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [0 .7 .9])
    plot([0 0], [-3 3], 'k-')
    plot([-3 3], [0 0], 'k-')
    plot([-3 3], [-3 3], 'k--')
    plot([-3 3], [3 -3], 'k--')
    axis([-3 3 -3 3])
    legend('Context 1, low contrast', 'Context 1, high contrast', 'Context 2, low contrast', 'Context 2, high contrast', 'location', 'NorthWest')
    title(strcat("Task design ", monkeyName))
    xlabel("Category bias (log odds)")
    ylabel("Direction bias (log odds)")
    
    subplot(2,4,3)
    plot(catLogOddsLC(dataSets,1), dirLogOddsLC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [1 .75 .75])
    hold on, box off, axis square
    plot(catLogOddsHC(dataSets,1), dirLogOddsHC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [1 0 0])
    plot(catLogOddsLC(dataSets,2), dirLogOddsLC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [.5 .9 1])
    plot(catLogOddsHC(dataSets,2), dirLogOddsHC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [0 .7 .9])
    plot([0 0], [-3 3], 'k-')
    plot([-3 3], [0 0], 'k-')
    plot([-3 3], [-3 3], 'k--')
    plot([-3 3], [3 -3], 'k--')
    axis([-3 3 -3 3])
    legend('Context 1, low contrast', 'Context 1, high contrast', 'Context 2, low contrast', 'Context 2, high contrast', 'location', 'NorthWest')
    title(strcat("Choice behavior ambiguous stimuli ", monkeyName))
    xlabel("Category bias (log odds)")
    ylabel("Direction bias (log odds)")
    
    subplot(2,4,4)
    plot(taskDirLogOddsLC(dataSets,1), dirLogOddsLC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [1 .75 .75])
    hold on, box off, axis square
    plot(taskDirLogOddsHC(dataSets,1), dirLogOddsHC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [1 0 0])
    plot(taskCatLogOddsLC(dataSets,1), catLogOddsLC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [.5 .9 1])
    plot(taskCatLogOddsHC(dataSets,1), catLogOddsHC(dataSets,1), 'ko', 'markersize', 8, 'markerfacecolor', [0 .7 .9])
    plot(taskDirLogOddsLC(dataSets,2), dirLogOddsLC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [1 .75 .75])
    plot(taskDirLogOddsHC(dataSets,2), dirLogOddsHC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [1 0 0])
    plot(taskCatLogOddsLC(dataSets,2), catLogOddsLC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [.5 .9 1])
    plot(taskCatLogOddsHC(dataSets,2), catLogOddsHC(dataSets,2), 'ko', 'markersize', 8, 'markerfacecolor', [0 .7 .9])    
    plot([0 0], [-2 2], 'k-')
    plot([-2 2], [0 0], 'k-')
    plot([-2 2], [-2 2], 'k--')
    axis([-2 2 -2 2])
    legend('Direction bias, low contrast', 'Direction bias, high contrast', 'Category bias, low contrast', 'Category bias, high contrast', 'location', 'NorthWest')
    title(strcat("Relation between task and behavior under ambiguity ", monkeyName))
    xlabel("Task design (log odds)")
    ylabel("Choice behavior (log odds)")
    
    subplot(2,4,7)
    plot(log2(oriSensLC(dataSets)), min(3.75, oriBiasLC(dataSets)), 'kd', 'markersize', 14, 'markerfacecolor', [1 1 1])
    hold on, box off, axis square
    plot(log2(oriSensHC(dataSets)), min(3.75, oriBiasHC(dataSets)), 'wd', 'markersize', 14, 'markerfacecolor', [0 0 0])
    plot(log2([oriSensLC(dataSets); oriSensLC(dataSets)]), min(3.75, oriBiasLC_IQR(dataSets,:)'), 'k-')
    plot(log2([oriSensHC(dataSets); oriSensHC(dataSets)]), min(3.75, oriBiasHC_IQR(dataSets,:)'), 'k-')
    plot(log2(oriSensLC_IQR(dataSets,:)'), min(3.75, [oriBiasLC(dataSets); oriBiasLC(dataSets)]), 'k-')
    plot(log2(oriSensHC_IQR(dataSets,:)'), min(3.75, [oriBiasHC(dataSets); oriBiasHC(dataSets)]), 'k-')
    plot(log2([.5 2]), [0 0], 'k--')
    axis([log2(.35) log2(2) -1 4])
    legend(sprintf('Low contrast, r = % 0.2g, P = %0.2g', rSensBiasLC, pSensBiasLC), sprintf('High contrast, r = % 0.2g, P = %0.2g', rSensBiasHC, pSensBiasHC), 'location', 'NorthWest')
    title(sprintf('r = % 0.2g, P = %0.2g (Spearman)', rSensBiasAll, pSensBiasAll))
    xlabel("Orientation sensitivity (log_2 basis)")
    ylabel("Orientation bias (deg)")
    
    % Plot example psychometric functions
    cd(dataPath)
    load(fileNameList(plotExpPrior(iA)).name);
    cd(functionPath)
    
    for iC = 1:numel(S.general.contrastValues)
        
        % Computations for plotting
        nSamples    = 100;
        stimPlot    = linspace(min(S.general.oriValuesDeg), max(S.general.oriValuesDeg), nSamples);
        [~, predPF] = giveNLL(S.beh.model(1).params{iC}, stimPlot, zeros(2,nSamples), zeros(2,nSamples), 'M1');
        
        % Plot example dataset
        if iC == 1
            subplot(2,4,6)
        else
            subplot(2,4,5)
        end
        
        plot(stimPlot, predPF(1,:), '-', 'linewidth', 2, 'color', [1 0 0])
        hold on, box off, axis square
        plot(stimPlot, predPF(2,:), '-', 'linewidth', 2, 'color', [0 .7 .9])
        plot([-5 5], [.5 .5], 'k--')
        plot([0 0], [0 1], 'k--')
        axis([-5 5 0 1])
        xlabel("Stimulus orientation (deg)")
        ylabel("Proportion clockwise")
        
        if iC == 1
            title("Low stimulus contrast")
        else
            title("High stimulus contrast")
        end
        
        
        for iO = 1:numel(S.general.oriValuesDeg)
            plot(S.general.oriValuesDeg(iO), S.beh.propCW_prior{iC}(1,iO), 'ko', 'markerfacecolor', [1 0 0], 'markersize', round(S.beh.nTrials_prior{iC}(1,iO)/5)+1)
            plot(S.general.oriValuesDeg(iO), S.beh.propCW_prior{iC}(2,iO), 'ko', 'markerfacecolor', [0 .7 .9], 'markersize', round(S.beh.nTrials_prior{iC}(2,iO)/5)+1)
        end
        legend('Context 1', 'Context 2', 'location', 'NorthWest')
    end
    cd(thisPath)
    
end







