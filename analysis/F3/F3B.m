% => context colors
clrs = {[1,0,0],[0,0.75,1]};
% => monkey session ranges
rgs = {1:13,14:29};
% => monkeys
M = {'F','J'};

% ==> only plot low contrast cases for F3B
rc = 1;

% ==> draw figure
figure();  set(gcf,'Color','w'); set(gcf, 'Position',[675 553 963 409])

for m = 1:length(M)
    subplot(1,2,(length(M) + 1)-m); % => switch order to match manuscript panel order
    for xc = 1:length(cx)
        hold on; hold all;
        % ==> plot average over model dv fit curves
        plot(linspace(-790,-40,200), mean(vertcat(dvsm{rgs{m},xc,rc}),1), 'color', clrs{xc},  'linewidth',4)
        % ==> show scatterplot of raw dv values within the curve fit time window
        % ==> timeSac(16) == -790 and timeSac(31) == -40
        rvs = mean(vertcat(dvsr{rgs{m},xc,rc}),1);
        scatter(timeSac(16:31), rvs(16:31), 75, 'o', 'markerfacecolor', clrs{xc},'markeredgecolor','w');
    end
    % ==> time windows
    % ==> delimiters for early dv values
    plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
    plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')
    % ==> delimiters for late dv values
    plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
    plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
    xlabel('Time');
    ylabel('Signed DV average (vertical stimulus)');    
    xlim([-800,0]);
    ylim([-0.2,0.2]);    
    title(['Monkey ',M{m}]);
end