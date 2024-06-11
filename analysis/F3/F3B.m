
figure();  set(gcf,'Color','w'); set(gcf, 'Position',[675 553 963 409])
subplot(1,2,1);
hold on; hold all;
p1 = plot(linspace(-790,-40,200), mean(J_cw_lo,1), 'color', [1,0,0,0.35],  'linewidth',8);
p3 = plot(linspace(-790,-40,200), mean(J_ccw_lo,1), 'color', [0,0,1,0.35], 'linewidth',8);
% ==> timeSac(16) == -790 and timeSac(31) == -40
p4 = plot(timeSac(16:31), mean(J_cw_lo_raw_DVs(:,16:31),1), 'o', 'color', [1,0,0,0.35]);
p5 = plot(timeSac(16:31), mean(J_ccw_lo_raw_DVs(:,16:31),1), 'o', 'color', [0,0,1,0.35], 'linewidth',1);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')

plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
xlim([-840,-30])
title('Monkey J');
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');
% legend([p1,p3],'cw context, low contrast', ...
%                'ccw context, low contrast');

subplot(1,2,2);
hold on; hold all;
p1 = plot(linspace(-795,-45,200), mean(F_cw_lo,1), 'color', [1,0,0,0.35],  'linewidth',8);
p3 = plot(linspace(-795,-45,200), mean(F_ccw_lo,1), 'color', [0,0,1,0.35], 'linewidth',8);
% ==> timeSac(16) == -790 and timeSac(31) == -40
p4 = plot(timeSac(16:31), mean(F_cw_lo_raw_DVs(:,16:31),1), 'o', 'color', [1,0,0,0.35]);
p5 = plot(timeSac(16:31), mean(F_ccw_lo_raw_DVs(:,16:31),1), 'o', 'color', [0,0,1,0.35], 'linewidth',1);
% => time windows
plot(repmat(-800,[20,1]),linspace(-0.2,0.2,20),'k--')
plot(repmat(-600,[20,1]),linspace(-0.2,0.2,20),'k--')

plot(repmat(-500,[20,1]),linspace(-0.2,0.2,20),'k:')
plot(repmat(-300,[20,1]),linspace(-0.2,0.2,20),'k:')
xlim([-840,-30])
title('Monkey F');
xlabel('Time');
ylabel('Signed DV average (vertical stimulus)');
legend([p1,p3],'cw context, low contrast', ...
               'ccw context, low contrast');
                 