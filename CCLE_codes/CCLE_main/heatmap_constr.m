function hmap = heatmap_constr(L,N,ts_traj,smc_lsites_traj,smc_rsites_traj,avgn,H,L0,kl)
% H is number of base pairs per simulation lattice. H = bpnum/L;
% iter=floor(length(ts_traj)/intv);
% kl is the Kuhn length of the polymer, usually = 1
% L0 is the size of the box

load('DSBZ2.mat')
ts_traj = ts_traj - ts_traj(1); % re-shift time steps after head elimination

bp_per_nano = 50;

lattice_to_kl = H / (kl * bp_per_nano); % 50 is number of basepair per nanometer of chromatin polymer. (# kuhn length per lattice)

% Calculate the look-up table based on the box_size and 
% lookup = contact_prob_box(linspace(0.01,100,10000),L0,kl); % look up table is up to 1k KL (good for most cases)
lookup_table = contact_prob_sphere_refbc(linspace(0.01,1000,100000),DSBZ2,L0,kl); % sphere: L0 is the sphere radius
% vq = interp1(linspace(0.1,1000,10000),lookup,linspace(0.01,1000,100000)); % bin is 0.01 Kuhn length
% vq = lookup_table;


Cor_NP=zeros(L,L);
Cor_NPa=Cor_NP;
for i = 1:avgn
    tm = ts_traj(end)/avgn*i;
    mid=find(ts_traj<=tm,1,'last');
    l_foot = smc_lsites_traj(mid,:);
    r_foot = smc_rsites_traj(mid,:);   
    for j=1:L-1
        for k=j+1:L
%             PT1=j;
%             PT2=k;            
%             l_foot = smc_lsites_traj(mid,:);
%             r_foot = smc_rsites_traj(mid,:);            
            Cor_NP(j,k) = SpatialCor_TwoPoints(l_foot,r_foot,j,k,N); % in unit of lattice site
        end
    end
    Cor_M = Cor_NP+Cor_NP';
%     Cor_M(isinf(Cor_M)) = 1e3;
    Cor_M(Cor_M <= 1) = 1; % if two lattices are les than 1 lattice (e.g. loop base), they are separated by 1 lattice
    Cor_M = Cor_M * lattice_to_kl; % convert unit from lattice number to Kuhn length

    % Input to look-up table
    lookup_input = round(Cor_M .* 100);
    Cor_M2 = lookup_table(lookup_input);
    
    if isnan(Cor_M2)
        disp(i)
    end
    Cor_NPa = Cor_NPa + Cor_M2;
    disp(i)
end
Cor_NPa = Cor_NPa/avgn;
Cor_NPa(logical(eye(L)))=0;
hmap=Cor_NPa;
end
