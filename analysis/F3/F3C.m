% ==> histograms for F3C

Fvs = cell(1,2);
Jvs = cell(1,2);

% ==> store model fits (within appropriate time window)
for xc = 1:length(cx)
    % ==> all dv model fit values
    Fv = vertcat(dvsm{1:13,xc,:});
    % ==> Monkey F all values for vertical stimulus - across hi and lo contrasts       
    Fvs{xc} = mean(Fv(:, il2:iu2),2); 
 
    % ==> all dv model fit values
    Jv = vertcat(dvsm{14:29,xc,:});   
    % ==> Monkey J all values for vertical stimulus - across hi and lo contrasts    
    Jvs{xc} = mean(Jv(:, il2:iu2),2); 
end

% ==> all data 
vs = {Fvs;Jvs};

% => contexts
cxp = {'ccw','cw'};
% => Monkeys
M = {'F','J'};
% => Monkey colors
clrs = {[1,0.75,0],[0.15,0.75,0.5]};

% ==> edges
edgs = -1.5:0.1:1.5;

fg = figure(); set(fg,'color','white'); set(fg, 'Position',[675 421 660 541])

sbp = 1;
for m = 1:length(M)
    for xc = 1:length(cx)
        subplot(2,2,sbp);
        h = histogram(vs{m}{xc},edgs,'facecolor',clrs{m}, 'edgecolor','w'); hold on; hold all;
        plot(zeros(1,max(h.Values)),1:max(h.Values),'r--','linewidth',2)
        title(['Animal ', M{m}, ': ', cxp{xc}, ' prior context'])
        xlabel('Signed DV initial offset')
        ylabel('Frequency')
        drawnow;
        sbp = sbp + 1;
        
        % ==> output wilcoxon signed test results
        [p,~] = signrank(vs{m}{xc});
        fprintf('Wilcoxon test for animal %s, %s cases, p = %d...\n',M{m}, cxp{xc}, p);
    end
    % ==> ranksum test
    [p,~] = ranksum(vs{m}{2}, vs{m}{1});
    fprintf('ranksum test for animal %s, p = %d...\n\n',M{m}, p);
end