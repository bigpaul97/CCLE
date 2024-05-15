classdef LEFSystem < handle
    properties
        L
        N
        time
        vels
        lifespans
        rebinding_times
        perms
        lattice
        locs
    end
    methods
        function obj = LEFSystem(L,N,vels,lifespans,rebinding_times,init_locs,BE_perms)
            
            % the occupancy of each site of the system. If -1, the site is unoccupied.
            obj.L = L;
            obj.N = N;
            obj.lattice = -1 * ones(L,1);
            obj.locs = -1 * ones(2*N,1);
            obj.vels=vels;
            obj.lifespans=lifespans;
            obj.rebinding_times = rebinding_times;                        
            
%             if isempty(BE)
%                 if isempty(perms)
%                     obj.perms(:,1) = ones(L,1);
%                     obj.perms(:,2) = ones(L,1);
%                 else
%                     obj.perms(:,1) = perms;
%                     obj.perms(:,2) = perms;
%                 end
%                 obj.perms(end,1)=0;
%                 obj.perms(1,end)=0;
%             else
%                 if isempty(perms)
%                     obj.perms(:,1) = ones(L,1); % non-CTCF region has perm = 1
%                     obj.perms(:,2) = ones(L,1);
%                 else
%                     obj.perms(:,1) = perms; % if a universal perms for non-CTCF region is assigned
%                     obj.perms(:,2) = perms;
%                 end
                % inherit perms to the class
            obj.perms(:,1) = BE_perms(:,1); % right barrier
            obj.perms(:,2) = BE_perms(:,2); % left barrier
%             end

%             obj.perms(1)=0;
%             obj.perms(end)=0;
%             
            % Initialize non-random loops
            for i = 1:2*N
                % Check if the loop is preinitialized.
                if init_locs(i) < 0
                    continue
                end
                % Populate a site.
                obj.locs(i) = init_locs(i);
                obj.lattice(obj.locs(i)) = i;
            end
        end
        
        function r = make_step(obj,leg_idx,direction)
            % The variable `direction` can only take values +1 or -1.
            new_pos = obj.locs(leg_idx) + direction;

            % make boundaries correct
            if new_pos == 0 % left bdy
                new_pos = obj.L;
            elseif new_pos > obj.L % right bdy
                new_pos = mod(new_pos, obj.L);
            end
            
            r = obj.move_leg(leg_idx, new_pos);
        end
        
        function r = move_leg(obj,leg_idx,new_pos)
            if new_pos > 0 && obj.lattice(new_pos) > 0
                r = 0;
                return
            end
            prev_pos = obj.locs(leg_idx);
            obj.locs(leg_idx) = new_pos;
            
            if prev_pos > 0
                if obj.lattice(prev_pos) <= 0
                    r = 0;
                    return
                end
                obj.lattice(prev_pos) = -1;
            end
            if new_pos > 0
                obj.lattice(new_pos) = leg_idx;
            end
            r = 1;
        end
        
        function okay = check_system(obj)
            okay = 1;
            for i = 1:obj.N
                if obj.locs(i) == obj.locs(i+obj.N)
                    disp(strcat('loop ' , num2str(i), 'has both legs at ', num2str(obj.locs(i))))
                    okay = 0;
                end
                if obj.locs(i) > obj.L
                    disp(strcat('leg ', num2str(i), 'is located outside of the system: ', num2str(obj.locs(i))))
                    okay = 0;
                end
                if obj.locs(i+obj.N) > obj.L
                    disp(strcat('leg ', num2str(i+obj.N), 'is located outside of the system: ', num2str(obj.locs(i+obj.N))))
                    okay = 0;
                end
                if (obj.locs(i) <= 0 && obj.locs(i+obj.N) > 0 ) || (obj.locs(i) > 0 && obj.locs(i+obj.N) <= 0 )
                    disp(strcat('the legs of the loop', num2str(i), 'are inconsistent: ', num2str((obj.locs(i))), num2str(obj.locs(i+obj.N))))
                    okay = 0;
                end
            end
        end
    end
end