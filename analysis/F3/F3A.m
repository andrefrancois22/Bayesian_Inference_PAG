
% ==> choose a good example (iS == 23 in the manuscript)
iS = 23;
% ==> load data for that one example
fprintf('loading data for session %d... \n',iS);
load([dataPath,'/dataSet_',num2str(iS),'.mat']);
fprintf('finished loading data! \n')

% => context colors
clrs = {[1,0,0],[0,0.75,1]};
% => monkey session ranges
rgs = {1:13,14:29};
% => monkeys
M = {'F','J'};

figure(); set(gcf,'Color','w'); set(gcf, 'Position',[675 363 602 599])

m = 2;
rc = 1;
for xc = 1:length(cx)
    hold on; hold all;
    % ==> plot average over model dv fit curves
    plot(linspace(-790,-40,200), mean(vertcat(dvsm{iS,xc,rc}),1), 'color', clrs{xc},  'linewidth',4)
    % ==> show scatterplot of raw dv values within the curve fit time window
    % ==> timeSac(16) == -790 and timeSac(31) == -40
    rvs = mean(vertcat(dvsr{iS,xc,rc}),1);
    scatter(timeSac(16:31), rvs(16:31), 75, 'o', 'markerfacecolor', clrs{xc},'markeredgecolor','w');   
end
% ==> time windows
% ==> delimiters for early dv values
plot(repmat(-800,[20,1]),linspace(-1,1,20),'k--')
plot(repmat(-600,[20,1]),linspace(-1,1,20),'k--')
% ==> delimiters for late dv values
plot(repmat(-500,[20,1]),linspace(-1,1,20),'k:')
plot(repmat(-300,[20,1]),linspace(-1,1,20),'k:')
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');    
xlim([-800,0]);
ylim([-1,1]); 
title([S.general.monkey,S.general.expDate])
ylim([-1,1]);
drawnow; 
