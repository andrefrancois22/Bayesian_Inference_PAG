
% => monkey colors
clrs = {[1,0.75,0],[0.15,0.75,0.5]};
% => monkey session ranges

rgs = {1:13,14:29};
% => monkeys
M = {'F','J'};
% => plot y axis ranges
xyl = {[-0.35,0.35],[-0.6,0.6]};

% => draw figure
figure(); set(gcf,'Color','w'); set(gcf,'Position',[178 358 1325 604]);

for m = 1:length(M)    
    subplot(1,2,m); hold on; hold all;
    for xc = 1:length(cx)
        for rc = 1:length(cr)
            % ==> scatter early and late dv values
            v = vertcat(dvw{rgs{m},xc,rc});
            scatter(v(:,1),v(:,2), 140, 'o','markerfacecolor', clrs{m}, 'markeredgecolor', 'w');
        end
    end
    % ==> lines
    plot(-0.6:0.1:0.6,-0.6:0.1:0.6,'k--')
    plot(-0.6:0.1:0.6,zeros([1,length(-0.6:0.1:0.6)]),'k--')
    plot(zeros([1,length(-0.6:0.1:0.6)]),-0.6:0.1:0.6,'k--')
    axis equal; axis square;
    xlabel('Average (-800ms to -600ms) window');
    ylabel('Average (-500ms to -300ms) window');
    title(['Monkey ',M{m}]);   
    xlim(xyl{m}); 
    ylim(xyl{m});      
end

% ==> eigenvectors for each monkey
Vs  = nan(2,2);
MUs = nan(2,2);
THs = nan(1,2);

% ==> plot first eigenvectors
for m = 1:length(M)  
    % ==> covariance matrix
    cv = cov(vertcat(dvw{rgs{m},:,:}));

    % ==> eigendecomposition
    [V,D] = eig(cv);
    % ==> get index of max eigenvalue
    [~,I] = max(D);
    % ==> mean of the distribution
    mu = mean(cv)';
    % ==> angle of first eigenvector
    [TH,~] = cart2pol(V(1,I(2)),V(2,I(2)));  
    THs(m) = TH;
    % ==> store Vs
    Vs(:,m) = V(:,I(2));
    % ==> store means
    MUs(:,m) = mu;
    
    subplot(1,2,m); hold on; hold all;
    % ==> plot eigenvector
    plot([mu(1),V(1,I(2))],[mu(2),V(2,I(2))],'r-','linewidth', 3) 
    title(['monkey ',M{m}, ' 1st PC slope (deg): ', num2str(rad2deg(THs(m)))]);
end
