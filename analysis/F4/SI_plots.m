% ==> draw Psychometric functions and proportions over PFs

% => number of sessions
Sn = 29;

for rc = 1:length(cr) % ==> show for lo or high contrast cases
    
    % ==> draw a new figure for each contrast condition
    figure(); set(gcf,'color','white'); set(gcf,'Position',[86 282 921 668]);
   
    % ==> what is the session
    for iS = 1:Sn
        
        % ==> orientations
        or = oris{iS};

        % ==> low dynamic range (and low contrast)
        % ==> count trials in ccw condition (==> <choice>_<context>_<contrast>_<dynamic range>)
        tn_ccw_lo = [cnts{2,1,rc,1,:}] + [cnts{1,1,rc,1,:}]; 
        % ==> count trials in cw condition
        tn_cw_lo =  [cnts{2,2,rc,1,:}] + [cnts{1,2,rc,1,:}];   
        % ==> high dynamic range (and low contrast)
        % ==> count trials in ccw condition (==> <choice>_<context>_<contrast>_<dynamic range>)
        tn_ccw_hi = [cnts{2,1,rc,2,:}] + [cnts{1,1,rc,2,:}];
        % ==> count trials in cw condition
        tn_cw_hi =  [cnts{2,2,rc,2,:}] + [cnts{1,2,rc,2,:}];   

        % ==> trial proportions for plot (just count, for each stimulus type, number of trials over total session trials)
        tn_ccw_lo_p = tn_ccw_lo ./ sum(tn_ccw_lo + tn_cw_lo); 
        tn_cw_lo_p =  tn_cw_lo  ./ sum(tn_ccw_lo + tn_cw_lo); 
        % ==> high dynamic range
        tn_ccw_hi_p = tn_ccw_hi ./ sum(tn_ccw_hi + tn_cw_hi);
        tn_cw_hi_p =  tn_cw_hi  ./ sum(tn_ccw_hi + tn_cw_hi);  

        % ==> figure specs
        line = {'.-',':'}; clr = {'b','r'; 'b','m'};
        gcf;    
        ax = subplot(5,6,iS); hold on; hold all;
        title(Ss{iS}); 
        hold on; hold all;  
        % => plot
        for xc = 1:length(cx)
            for sp = 1:2
                plot(or,PF{iS,rc,sp}(xc,:),[clr{rc,xc},line{sp}],'linewidth', 1)    
            end
        end
        plot(zeros(length(or),1),linspace(0,1,length(or)),'k--','linewidth', 0.5) 
        plot(or,ones(length(or),1)/2,'k--','linewidth', 0.5) 
        ax.XTick = or;
        ax.XTickLabel = or;
        xlabel('orientation')
        ylabel('prop cw')  

        % ==> for each orientation plot proportions
        for o = 1:length(or)
            scatter(or(o), CPs{iS,rc,1}(2,o),  max([1,round(tn_ccw_lo_p(o)*300)]), 'o','filled', 'markerfacecolor', [1,0,0], 'markeredgecolor', clr{rc,2})
            scatter(or(o), CPs{iS,rc,1}(1,o),  max([1,round(tn_cw_lo_p(o)*300)]),  'o','filled', 'markerfacecolor', [0,0,1], 'markeredgecolor', clr{rc,1})
            % ==> trial proportions for high dynamic range
            scatter(or(o), CPs{iS,rc,2}(2,o),  max([1,round(tn_ccw_hi_p(o)*300)]), 'o','filled', 'markerfacecolor', [1,0.5,0.5], 'markeredgecolor', clr{rc,2})
            scatter(or(o), CPs{iS,rc,2}(1,o),  max([1,round(tn_cw_hi_p(o)*300)]),  'o','filled', 'markerfacecolor', [0.5,0.5,1], 'markeredgecolor', clr{rc,1})
        end
        axis square;
        drawnow;  
    end
end