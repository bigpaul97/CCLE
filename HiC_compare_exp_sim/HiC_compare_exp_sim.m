function [hic_compare_map,ps_exp,mpr,ps_mpr,pearson_coeff] = ...
    HiC_compare_exp_sim(filename,exphmap,avenum,L)
% Generate Hi-C map comparison between experiment and simulation, P(s) plot, ratio map, and loop occupancy plot

% Inputs
% filename: the MATLAB filename with the simulated Hi-C map (without number indices at the end)
% exphmap: experimental Hi-C map for comparison
% avenum: range of simulation file indices to average, e.g. 1:200
% L: number of lattice sites in the simulation

% Outputs
% hic_compare_map: the matrix of averaged simulated map (lower triangle) versus experimental map (upper triangle)
% ps_exp: P(s) of experimental Hi-C
% mpr: Mean Pairwise Ratio (MPR)
% ps_mpr: Mean Pairwise Ratio for P(s) (MPR_{P(s)})
% pearson_coeff: Pearson Correlation Coefficient (PCC)


%% Initialization

exphmap_dim = size(exphmap,1); 

% range of the experimental Hi-C to include to calculate the score. The plotting always uses the whole range
bdy_range = 26:275; % pombe:11:110 (10k), 51:550 (2k); meiotic scere: 26:275; mitotic scere: 301:500

% color scale of the Hi-C map comparison. Can be changed to improve visualization
cmap_low_limit = -2.699; % pombe: -2.5229;-2.58(2k), scere: -2.6990
cmap_high_limit = -0.699; % pombe: -0.5229;-0.58(2k), scere: -0.6990

% Implement Gaussian scaling or not
gauss_YN = 1; % 1: yes; 0: no
gauss_sigma = 0.9e5; % only factor matters is the SD of the gaussian

% plot the loop occupancy or not (relatively slow, if loop occupancy info is not needed, spare it)
loop_config_YN = 1; % 1: yes; 0: no. 

% specify the diagonal whose P(s) should match that of experiment
diag_scaling_index = 2; % pombe: 2; meiotic scere: 2; mitotic scere: 3

% Number of diagonals of the experimental Hi-C to include to calculate metric scores
Ndiagonal = 32; % pombe: 12 (10k-resln); 60 (2k-resln). meiotic scere: 32 (2k-resln). mitotic scere: 90 (500-resln).

% number of points to exclude for P(s).
tail_crop = 240; % pombe: 70 (10k-resln); 350 (2k-resln). meiotic scere: 240 (2k-resln). mitotic scere: 1000 (500-resln).


%% Average and bin simulation Hi-C map 

simhmap = aveHiC(avenum,filename); % change this when compare exp vs exp

simhmap_dim = size(simhmap,1);
binf = simhmap_dim/exphmap_dim;

% Bin simulation hic to the exphic resolution
dg = -diag(simhmap) + 1;
simhmap = simhmap + diag(dg,0);
fun = @(block_struct) mean(block_struct.data(:));
B = blockproc(simhmap, [binf binf], fun); % simhmap is binned to the same resolution as exphic
B = B./max(reshape(B,[],1));

% set hmap_sim diagonal and first diagonals to zeros
hmap_sim = B - diag(diag(B)) - diag(diag(B,1),1) - diag(diag(B,-1),-1); 


%% Gaussian factor

if gauss_YN

    load(append(filename,'_1','.mat'),'bpnum')
    
    gauss_array_dim = exphmap_dim - 1; % there are (exphmap_dim - 1) diagonals
    
    gauss_sep = bpnum/(gauss_array_dim + 1);
    gauss_x = gauss_sep:gauss_sep:gauss_array_dim * gauss_sep;
    gauss_value = exp(-(gauss_x-0).^2/2/(gauss_sigma)^2);
    
    % Scale each diagonal by the corresponding factor from Gaussian distr
    hmap_sim_gauss_scaled = zeros(exphmap_dim,exphmap_dim);
    for i = 1:gauss_array_dim
        hmap_sim_gauss_scaled = hmap_sim_gauss_scaled + gauss_value(i) * diag(diag(hmap_sim,i),i);
    end
    
    hmap_sim = hmap_sim_gauss_scaled + hmap_sim_gauss_scaled'; % The new hmap_sim is gauss scaled

end


%% Plot the P(s) curve
load(append(filename,'_1','.mat'),'H')
ps_bin = H * binf; % # of basepair per exphmap lattice

% set exphmap diagonal and first diagonals to zeros
exphmap = exphmap - diag(diag(exphmap)) - diag(diag(exphmap,1),1) - diag(diag(exphmap,-1),-1);

ps_exp = zeros(exphmap_dim-1,1);
ps_sim = zeros(exphmap_dim-1,1);

% P(s) for non-zero pixels
for i = 1:exphmap_dim-1
    diag_exp = diag(exphmap,i);
    diag_sim = diag(hmap_sim,i);

    ps_exp(i) = mean(diag_exp);
    ps_sim(i) = mean(diag_sim);
end

% make sim to match exp initial values (for both hic map and p(s) curve)
sim_corrf = ps_exp(diag_scaling_index)/ps_sim(diag_scaling_index); % change accordingly
hmap_sim = hmap_sim .* sim_corrf;
ps_sim = ps_sim .* sim_corrf;

% p(s) mean pairwise ratio (exclude the first and second order gap diagonals)
ps_mpr = mean(exp(abs(log(ps_exp(3:end-tail_crop)) - log(ps_sim(3:end-tail_crop)))));

% Plot
figure(), hold on

scatter(ps_bin*2:ps_bin:ps_bin*(exphmap_dim-1),ps_exp(2:end),'o','linewidth',1,'MarkerEdgeColor','#0072bd') % 0072BD
plot(ps_bin*2:ps_bin:ps_bin*(exphmap_dim-1),ps_sim(2:end),'linewidth',2,'DisplayName','Sim','color',"#d95319") % D(5319
set(gca,'XScale','log','YScale','log','fontsize',20)
set(gca,'LineWidth',1.2)
xlabel('Genomic separation s (base pair)')
ylabel('Mean contact probability P(s)')

box on
hold off

%% Plot Hi-C comparison

% plot the Hi-C comparison
hic_compare_map = trianglize(exphmap,hmap_sim,exphmap_dim); % put exp upper-right, sim lower-left

figure()
imagesc(log10(hic_compare_map))
colormap(flipud(hot(256)));
hcb1 = colorbar;
set(hcb1,'ticklength',0);
clim([cmap_low_limit cmap_high_limit]);
pbaspect([1 1 1]);
set(gca,'TickLength',[0 0])


%% Plot difference map

% Exclude hmap_sim pixels where exphmap data is missing
hmap_sim(exphmap == 0) = 0;

exphmap_log = log10(exphmap);
hmap_sim_log = log10(hmap_sim);
exphmap_log(isinf(exphmap_log)) = 0;
hmap_sim_log(isinf(hmap_sim_log)) = 0;

diff_map = abs(hmap_sim_log - exphmap_log);

figure()
imagesc(diff_map)
colormap(flipud(gray(256)))
hcb2 = colorbar;
set(hcb2,'ticklength',0);
clim([0 0.6021]);
pbaspect([1 1 1]);
set(gca,'TickLength',[0 0])


%% getting the Hi-C score (the diagonal, first diagonals, edges and corners are excluded)

hic_score_matrix = diff_map;

% Only keep data upto number of diagonals given by Ndiagonal
y_size = size(hic_score_matrix,1);
x_size = size(hic_score_matrix,2);
hic_score_matrix_nocorner = zeros(y_size,x_size);

for i = -Ndiagonal:Ndiagonal
    hic_score_matrix_nocorner = hic_score_matrix_nocorner + diag(diag(hic_score_matrix,i),i);
end

% Get boundaries of the maps removed
y_range_elim = bdy_range;
x_range_elim = bdy_range;

hic_score_matrix_elim = hic_score_matrix_nocorner(y_range_elim,x_range_elim);


%% Calculate loop density score

if loop_config_YN
    stablePt = 1.5e4; % only steady-state loop configurations count
    [simoccup,~] = CS_vs_SiteOccup(filename,avenum,stablePt,L);
end

%% Pearson score between exphic and simhic (log-value)

% exphic
exphic_ps_scaled = zeros(y_size,x_size);
for i = 2:Ndiagonal
    dd = diag(exphmap_log,i);
    dd = 10.^dd;
    dd(dd == 1) = 0;
    dd_norm = dd./mean(dd);
    exphic_ps_scaled = exphic_ps_scaled + diag(dd_norm,i);
end

exphic_ps_scaled = exphic_ps_scaled + exphic_ps_scaled';
exphic_ps_scaled = exphic_ps_scaled(y_range_elim,x_range_elim);

% set the zero-far-diagonal entries to the average of the close-diagonal
avg_exp = mean(exphic_ps_scaled(exphic_ps_scaled ~= 0));
exphic_ps_scaled(exphic_ps_scaled == 0) = avg_exp;

% simhic
simhic_ps_scaled = zeros(y_size,x_size);
for i = 2:Ndiagonal
    dd = diag(hmap_sim_log,i);
    dd = 10.^dd;
    dd(dd == 1) = 0;
    dd_norm = dd./mean(dd);
    simhic_ps_scaled = simhic_ps_scaled + diag(dd_norm,i);
end

simhic_ps_scaled = simhic_ps_scaled + simhic_ps_scaled';
simhic_ps_scaled = simhic_ps_scaled(y_range_elim,x_range_elim);

avg_sim = mean(simhic_ps_scaled(simhic_ps_scaled ~= 0));
simhic_ps_scaled(simhic_ps_scaled == 0) = avg_sim;


%% Calculate combined score

hic_score_matrix_elim(hic_score_matrix_elim==0) = -Inf;
hic_exponentiate = 10.^hic_score_matrix_elim;

mpr = mean(hic_exponentiate(hic_exponentiate ~= 0));
pearson_coeff = corr2(exphic_ps_scaled,simhic_ps_scaled);



end