% ==> plot markers and colors
m = {'r.-','r+:','rh:','k^-'; 'bo-','bx:','bh:','kh:'};
c = {'r','b'};

% ==> accuracy plot
fig = figure(1); grid minor;
subplot(2,3,I); hold on; hold all; 
fig.Position = [249 177 1580 773]; 
set(gcf,'color','white');
plot(sigs,mle_v(I,:),'k^-', 'linewidth',1); hold on; hold all;
plot(sigs,map_v(I,:),m{I,1},'linewidth',1); hold on; hold all;
title(['Accuracy (prior context = ', ctx, ')'])
ylabel(['Overall accuracy (context = ', ctx, ')']);
xlabel('Sensory uncertainty \sigma');
legend('MLE','MAP',[ctx,' context ','posterior density > 0.5'], [ctx, ' context ', 'posterior sample'], 'Location','SouthWest')

subplot(2,3,3); grid minor;
plot(sigs,mean(mle_v,1),'k^-', 'linewidth',1); hold on; hold all;
plot(sigs,mean(map_v,1),m{1,1},'linewidth',1); hold on; hold all;
title('Overall Accuracy');
ylabel('Overall accuracy');
xlabel('Sensory uncertainty \sigma');
legend('MLE','MAP', 'Location','NorthEast')

idx = 1;
for i = [1,length(sigs)]
    sig = sigs(i);
    % ==> plot 
    d = (1:n_o) - vidx;
    subplot(2,3,3+idx); grid minor; hold on; hold all; 
    plot(d,squeeze(pr_mle(I,i,:))',m{I,4}, 'linewidth',0.5);
    plot(d,squeeze(pr_map(I,i,:))',m{I,1}, 'linewidth',0.5);      
    title(['Proportion "cw", \sigma = ', num2str(sig)]);
    xlabel('Orientations \theta');
    ylabel('proportion "cw" response');
    drawnow;
    % => subplot index
    idx = idx + 1;
end
subplot(2,3,5)
plot(d,predPF_map(1,:),'m-')
plot(d,predPF_map(2,:),'m-')
plot(d,predPF_mle(1,:),'g-')
plot(d,predPF_mle(2,:),'g-')
% => legend
legend('MLE cw context', ...
       'MAP cw context', ...
       'MAP cw PF fit', ...
       'MAP ccw PF fit', ...
       'MLE cw PF fit', ...
       'MLE ccw PF fit', ...
       'MLE ccw context', ...
       'MAP ccw context', ...
       'Location', 'SouthEast', 'Fontsize', 7);

subplot(2,3,6); grid minor;
hold on; hold all;
plot(sigs, bsv_map,'r.-');  
plot(sigs, bsv_mle,'k^-')
title('Bias and perceptual uncertainty')
ylabel('PF Bias');
xlabel('Sensory uncertainty \sigma');
legend('MAP','MLE', 'Location','NorthWest')
