function [siteOccup,y1axis] = CS_vs_SiteOccup(fileName,aveNum,stablePt,L)
% Calculate site occupation by fraction of time each locus is occupied by a leg
% The answer should be approximately the ChIP-Seq data for consistency
% (Note: time pops when new action is accomplished, so t_2-t_1 should be
% the duration of action 1)

siteOccup = zeros(L,1);

for i = aveNum
    load(append(fileName,'_',num2str(i),'.mat'),'l_sites_traj','r_sites_traj','ts_traj','CTCF','N','H')
    normalize = ts_traj(end) - ts_traj(stablePt); % normalize by the simulation time
    for j = stablePt:length(ts_traj)-1
        timestep = ts_traj(j+1) - ts_traj(j);
        a = l_sites_traj(j,:);
        b = r_sites_traj(j,:);
        siteOccup(a(a>0)) = siteOccup(a(a>0)) + timestep/normalize;
        siteOccup(b(b>0)) = siteOccup(b(b>0)) + timestep/normalize;
    end
    disp(i)
end

% take average over sample and over # of SMC legs
% siteOccup = siteOccup/aveNum/N/2;

% the head and tails are artifacts, set their occupancy to an average value
% Then renormalize total measure to 1.
% aa = mean(siteOccup(2:end-1));
% siteOccup(1) = aa;
% siteOccup(end) = aa;

% make the data format consistent
% if size(CTCF,2) == 1
%     CTCF(:,2) = CTCF(:,1);
% end

siteOccup = siteOccup/sum(siteOccup)*N*2; % normalized by number of LEF anchors

%% plot
y1axis = CTCF(:,2)/sum(CTCF(:,2))*N*2*(1000/H); % normalized by number of LEF anchors and resolution
y2axis = siteOccup*(1000/H); % normalized by resolution
xaxis = 1:L;

figure(), hold on
stairs(xaxis,y1axis,'linewidth',2,'DisplayName','Cohesin ChIP-seq data')
stairs(xaxis,y2axis,'linewidth',2,'DisplayName','Loop Occup')

xlabel('Lattice site','fontsize',32)
ylabel('Mean occupancy probability (1/kb)','fontsize',32)
% title('cohesin chip-seq vs simulated loop occupancy','fontsize',24)

y_limit = max(max(y1axis),max(y2axis))*1.2;

ylim([0 y_limit])

legend off
set(gca,'FontSize',32)
set(gca,'LineWidth',2)
set(gcf,'Position',[1,1,1321,737])
box on
hold off

end

