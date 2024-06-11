% ==> histograms for F3C
close all; clc;

% ==> all CW offset values for vertical stimulus - across hi and lo contrasts
cwsF  = [vcwlosF; vcwhisF];   
% ==> initial offset
cwsFi  = cwsF(:,1);
% ==> all CCW offset values for vertical stimulus - across hi and lo contrasts
ccwsF = [vccwlosF; vccwhisF]; 
% ==> initial offset
ccwsFi = ccwsF(:,1);

% ==> all CW offset values for vertical stimulus - across hi and lo contrasts
cwsJ  = [vcwlosJ; vcwhisJ];   
% ==> initial offset
cwsJi  = cwsJ(:,1);
% ==> all CCW offset values for vertical stimulus - across hi and lo contrasts
ccwsJ = [vccwlosJ; vccwhisJ]; 
% ==> initial offset
ccwsJi = ccwsJ(:,1);

% ==>
cwsF_mu  = mean(cwsF(:,il2:iu2),2);
ccwsF_mu = mean(ccwsF(:,il2:iu2),2);

cwsJ_mu  = mean(cwsJ(:,il2:iu2),2);
ccwsJ_mu = mean(ccwsJ(:,il2:iu2),2);

% ==> edges
edgs = -1.5:0.1:1.5;

fg = figure(); set(fg,'color','white'); set(fg, 'Position',[675 421 660 541])
subplot(2,2,1);
hFcw = histogram(cwsF_mu,edgs,'facecolor',[0.15,0.75,0.5]); hold on; hold all;
plot(zeros(1,max(hFcw.Values)),1:max(hFcw.Values),'r--','linewidth',2)
title('Animal F: cw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;
subplot(2,2,2);
hFccw = histogram(ccwsF_mu,edgs,'facecolor',[0.15,0.75,0.5]); hold on; hold all;
plot(zeros(1,max(hFccw.Values)),1:max(hFccw.Values),'r--','linewidth',2)
title('Animal F: ccw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;
 
subplot(2,2,3);
hJcw = histogram(cwsJ_mu,edgs,'facecolor',[1,0.5,0]); hold on; hold all;
plot(zeros(1,max(hJcw.Values)),1:max(hJcw.Values),'r--','linewidth',2)
title('Animal J: cw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;
subplot(2,2,4);
hJccw = histogram(ccwsJ_mu,edgs,'facecolor',[1,0.5,0]); hold on; hold all;
plot(zeros(1,max(hJccw.Values)),1:max(hJccw.Values),'r--','linewidth',2)
title('Animal J: ccw prior context')
xlabel('Signed DV initial offset')
ylabel('Frequency')
drawnow;

% ==> 
% ==> signrank wilcoxon test
[p_cwJ,~]  = signrank(cwsJ_mu);
[p_ccwJ,~] = signrank(ccwsJ_mu);

[p_cwF,~]  = signrank(cwsF_mu);
[p_ccwF,~] = signrank(ccwsF_mu);

fprintf('Wilcoxon test p values evaluating initial offsets in the vertical stimulus condition...\n')
fprintf('Wilcoxon test for animal F cw cases, p = %d...\n', p_cwF)
fprintf('Wilcoxon test for animal F ccw cases, p = %d...\n', p_ccwF)


fprintf('Wilcoxon test for animal J cw cases, p = %d...\n', p_cwJ)
fprintf('Wilcoxon test for animal J ccw cases, p = %d...\n', p_ccwJ)
 
% ==> ranksum
[PJ,~] = ranksum(cwsJ_mu,ccwsJ_mu)

[PF,~] = ranksum(cwsF_mu,ccwsF_mu)

% ==> restrict to the low contrast trials and repeat ranksum