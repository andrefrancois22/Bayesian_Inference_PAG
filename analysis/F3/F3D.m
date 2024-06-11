
figure(); set(gcf,'Color','w'); set(gcf,'Position',[178 358 1325 604]);
subplot(1,2,1);
hold on; hold all;
% => monkey F
s1 = scatter(Fmu_cw_lo(:,1), Fmu_cw_lo(:,2),  140, 'o','markerfacecolor', [0.25,0.85,0.6], 'markeredgecolor', 'r');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Fmu_ccw_lo(:,1),Fmu_ccw_lo(:,2), 140, 'o','markerfacecolor', [0.25,0.85,0.6], 'markeredgecolor', 'b');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s3 = scatter(Fmu_cw_hi(:,1), Fmu_cw_hi(:,2),  140, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'r');
s4 = scatter(Fmu_ccw_hi(:,1),Fmu_ccw_hi(:,2), 140, 'o','markerfacecolor', [0.15,0.75,0.5], 'markeredgecolor', 'b');
% ==> lines
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window');
ylabel('Average (-500ms to -300ms) window');
legend([s1,s2,s3,s4],'cw low contrast','ccw low contrast','cw high contrast','ccw high contrast','Location','SouthEast')
title('Monkey F');

subplot(1,2,2);
hold on; hold all;
% => monkey J
s1 = scatter(Jmu_cw_lo(:,1), Jmu_cw_lo(:,2),  140, 'o','markerfacecolor', [1,0.75,0], 'markeredgecolor', 'r');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s2 = scatter(Jmu_ccw_lo(:,1),Jmu_ccw_lo(:,2), 140, 'o','markerfacecolor', [1,0.75,0], 'markeredgecolor', 'b');%, 'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',1);
s3 = scatter(Jmu_cw_hi(:,1), Jmu_cw_hi(:,2),  140, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'r');
s4 = scatter(Jmu_ccw_hi(:,1),Jmu_ccw_hi(:,2), 140, 'o','markerfacecolor', [1,0.5,0], 'markeredgecolor', 'b');
% ==> lines
plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
axis equal; axis square;
xlabel('Average (-800ms to -600ms) window')
ylabel('Average (-500ms to -300ms) window')
legend([s1,s2,s3,s4],'cw low contrast','ccw low contrast','cw high contrast','ccw high contrast','Location','SouthEast')
title('Monkey J');