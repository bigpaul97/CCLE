function []=CCLE_main(tra)
% The CCLE simulation

%% Loading data
disp('Running CCLE simulation')

% Load the ChIP-seq data
% Note that the name of the example variables loaded is "CTCF" for some historical reason. 
% But here the variable represents the ChIP-seq data of the protein of interest
load('scere_cohesin_chipseq/chr13_rec8_240_840_500_shifted')

% random seed
rng(tra)


%% Adjustable parameters

% LEF simulation
bpnum = 6e5; % pombe: 1.2e6; yeast: 6e5
L = 1200; % Number of lattice sites (pombe: 1200; scere: 1200)
N = 40; % number of LEFs
vmulti = 5e-4; % The extrusion rate multiplier. This is the nominal rate, r, where r = rho*L/tau
N_SNAPSHOTS = 1e5; % Total number of LEF event simulated. 1e5 is enough, the longer, the smoother the map. 

% Hi-C map generation
intv = 100; % Number of averaged maps (should be 20000 maps in total)
steady = 1.5e4; % Estimated start of steady state (All LEF events before are discarded)
H = bpnum/L; % Number of basepair per simulationlattice
L0 = 500; % radius of the confinement sphere (nanometer)
kl = 160; % persistence length (nanometer)

%% Fixed parameters

R_EXTEND = 1; % This constant extrusion rate is fixed
R_SHRINK = 0; % Only outward steps allowed
R_OFF = 5e-4; % the dissociation rate is also fixed tau = 2000s 
T_MAX = 5e2; % a control parameter, don't change

% set zero or negative ChIP-seq signal to the minimum signal
minimum_signal = min(CTCF(CTCF(:,2) > 0,2)); % find the minimum non-zero signal
CTCF(CTCF(:,2) <= 0,2) = minimum_signal;
P_occup = CTCF(:,2) / sum(CTCF(:,2)) * 2 * N;

BE_perms = zeros(L,2); 
BE_perms(1:end-1,1) = 1 ./ (P_occup(1:end-1).*(1-P_occup(2:end))) * vmulti; % right barrier
BE_perms(2:end,2) = 1 ./ (P_occup(2:end).*(1-P_occup(1:end-1))) * vmulti; % left barrier

INIT_L_SITES = -1 * ones(N,1);
INIT_R_SITES = -1 * ones(N,1);
ACTIVATION_TIMES = zeros(N,1);
REBINDING_TIME = 0;
verbose = 1;


%% Initialization
VELS = zeros(4*N,1);
LIFESPANS = zeros(N,1);
REBINDING_TIMES = zeros(N,1);

for i = 1:N
    if length(R_EXTEND) > 1
        VELS(i) = R_EXTEND(i);
        VELS(i+3*N) = R_EXTEND(i);
    else
        VELS(i) = R_EXTEND;
        VELS(i+3*N) = R_EXTEND;
    end
    
    if length(R_SHRINK) > 1
        VELS(i+N) = R_SHRINK(i);
        VELS(i+2*N) = R_SHRINK(i);
    else
        VELS(i+N) = R_SHRINK;
        VELS(i+2*N) = R_SHRINK;
    end
    
    if length(R_OFF) > 1
        LIFESPANS(i) = 1/R_OFF(i);
    else
        LIFESPANS(i) = 1/R_OFF;
    end
    
    if length(REBINDING_TIME) > 1
        REBINDING_TIMES(i) = REBINDING_TIME(i);
    else
        REBINDING_TIMES(i) = REBINDING_TIME;
    end
end

INIT_LOCS = -1 * ones(2*N,1);

for i = 1:N
    INIT_LOCS(i) = INIT_L_SITES(i);
    INIT_LOCS(i+N) = INIT_R_SITES(i);
end

for i = 1:N
    if INIT_LOCS(i) ~= -1
        if INIT_LOCS(i+N) == -1
            warning('The initial positions of the right legs are chosen randomly, but those of the left legs are not.')
        end
        if ACTIVATION_TIMES(i) ~= 0
            warning('With defined INIT_L_SITES and INIT_R_SITES, ACTIVATION_TIMES must be 0.')
        end
    else
        if INIT_LOCS(i+N) ~= -1
            warning('The initial positions of the left legs are chosen randomly, but those of the right legs are not.')
        end
    end
end


LEFSYSTEM = LEFSystem(L, N, VELS, LIFESPANS, REBINDING_TIMES, INIT_LOCS,BE_perms);

l_sites_traj = zeros(N_SNAPSHOTS, N);
r_sites_traj = zeros(N_SNAPSHOTS, N);
ts_traj = zeros(N_SNAPSHOTS, 1);


prev_snapshot_t = 0;
snapshot_idx = 1;


evheap = Event_heap();

% Move LEFs onto the lattice at the corresponding activations times.
% If the positions were predefined, initialize the fall-off time in the
% standard way.
for i = 1:LEFSYSTEM.N
    % if the loop location is not predefined, activate it
    % at the predetermined time
    if INIT_LOCS(i) == -1
        evheap.add_event(i + 5 * LEFSYSTEM.N, ACTIVATION_TIMES(i));
        
        % otherwise, the loop is already placed on the lattice and we need to
        % regenerate all of its events and the motion of its neighbours
    else
        regenerate_all_loop_events(LEFSYSTEM, evheap, i)
        regenerate_neighbours(LEFSYSTEM, evheap, INIT_LOCS(i))
        regenerate_neighbours(LEFSYSTEM, evheap, INIT_LOCS(i+LEFSYSTEM.N))
    end
end

while snapshot_idx <= N_SNAPSHOTS

    LEFEvent = evheap.pop_event();
    LEFSYSTEM.time = LEFEvent.time;
    event_idx = LEFEvent.event_idx;
    
    LEFStatus = do_event(LEFSYSTEM, evheap, event_idx);
    
    if LEFStatus == 0
        disp('an assertion failed somewhere')
        return
    end
    
    if LEFSYSTEM.time > prev_snapshot_t + T_MAX / N_SNAPSHOTS
        prev_snapshot_t = LEFSYSTEM.time;
        l_sites_traj(snapshot_idx,1:N) = LEFSYSTEM.locs(1:N);
        r_sites_traj(snapshot_idx,1:N) = LEFSYSTEM.locs(N+1:end);
        ts_traj(snapshot_idx) = LEFSYSTEM.time;
        snapshot_idx = snapshot_idx + 1;
        if verbose && mod(snapshot_idx,1e4) == 0
            disp([snapshot_idx/N_SNAPSHOTS])
        end
    end
end

%% Hi-C map Generation

hmap = heatmap_constr(L,N,ts_traj(steady:end),l_sites_traj(steady:end,:),r_sites_traj(steady:end,:),intv,H,L0,kl);

%% Save data

savepath = '/home/ty268/project/loopEx/test_results/cohesin_LE_varvel/final_opt/chr2_200_1400_interphase_cohesin/';
fname = append('chiprescale40_LEF40_vmulti5e-4_kl160_a500','_',num2str(tra),'.mat');

save([savepath fname]);

end
