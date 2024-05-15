function sta = do_event(LEFSystem, evheap, event_idx)

%     Apply an event from a heap on the system and then regenerate it.
%     If the event is currently impossible (e.g. a step onto an occupied site),
%     it is not applied, however, no warning is raised.
%
%     Also, partially checks the system for consistency. Returns 0 if the system
%     is not consistent (a very bad sign), otherwise returns 1 if the event was a step
%     and 2 if the event was rebinding.
%
%     Possible events:
%     1 to 2N : a step to the left
%     2N+1 to 4N : a step to the right
%     4N+1 to 5N : passive unbinding
%     5N+1 to 6N : rebinding to a randomly chosen site

if event_idx <= 4 * LEFSystem.N
    % Take a step
    if event_idx <= 2 * LEFSystem.N
        leg_idx = event_idx;
        direction = -1;
    else
        leg_idx = event_idx - 2 * LEFSystem.N;
        direction = 1;
    end
    
    prev_pos = LEFSystem.locs(leg_idx);
    % check if the loop was attached to the chromatin
    sta = 1;
    if prev_pos > 0
        % make a step only if there is no boundary and the new position is unoccupied
        if LEFSystem.perms(prev_pos, (3-direction)/2) > 0
            if (prev_pos + direction) == 0
                new_position = LEFSystem.L;
            elseif (prev_pos + direction) > LEFSystem.L
                new_position = mod(prev_pos + direction, LEFSystem.L);
            elseif (prev_pos + direction) > 0 && (prev_pos + direction) <= LEFSystem.L
                new_position = prev_pos + direction;
            end

%             disp(new_position)

            if LEFSystem.lattice(new_position) < 0
                sta = sta * LEFSystem.make_step(leg_idx, direction);
                % regenerate events for the previous and the new neighbors
                regenerate_neighbours(LEFSystem, evheap, prev_pos); 
                regenerate_neighbours(LEFSystem, evheap, new_position); % CORRECT HERE!!
            end
        end
        
        % regenerate the performed event
        regenerate_event(LEFSystem, evheap, event_idx);
    end
    
elseif event_idx > 4 * LEFSystem.N && event_idx <= 5 * LEFSystem.N
    % unbinding 
    loop_idx = event_idx - 4 * LEFSystem.N;
    
    sta = 2;
    % check if the loop was attached to the chromatin
    if LEFSystem.locs(loop_idx) <= 0 || LEFSystem.locs(loop_idx+LEFSystem.N) <= 0
        sta = 0;
    end
    
    % save previous positions, but don't update neighbours until the loop
    % has moved
    % pos1 index could be larger than pos2 index
    prev_pos1 = LEFSystem.locs(loop_idx);
    prev_pos2 = LEFSystem.locs(loop_idx + LEFSystem.N);
    
    sta = sta * LEFSystem.move_leg(loop_idx, -1);
    sta = sta * LEFSystem.move_leg(loop_idx+LEFSystem.N, -1);
    
    % regenerate events for the loop itself and for its previous neighbours
    regenerate_all_loop_events(LEFSystem, evheap, loop_idx);
    
    % update the neighbours after the loop has moved
    regenerate_neighbours(LEFSystem, evheap, prev_pos1);
    regenerate_neighbours(LEFSystem, evheap, prev_pos2);
    
elseif event_idx > 5 * LEFSystem.N && event_idx <= 6 * LEFSystem.N
    loop_idx = event_idx - 5 * LEFSystem.N;
    
    sta = 2;
    % check if the loop was not attached to the chromatin
    if LEFSystem.locs(loop_idx) > 0 || LEFSystem.locs(loop_idx+LEFSystem.N) > 0
        sta = 0;
    end
    
    % find a new position for the LEM (a brute force method, can be
    % improved)
    while 1
        new_pos1 = randi([1 LEFSystem.L-1],1,1); % new LEF can't bind on the boundaries
        new_pos2 = new_pos1 + 1;
        if new_pos2 > LEFSystem.L
%             new_pos2 = mod(new_pos2,LEFSystem.L);
            continue
        end

        if LEFSystem.lattice(new_pos1) <= 0 && LEFSystem.lattice(new_pos2) <= 0
            break % could get into infinite loop if leg number is bigger than # lattice sites
        end
    end

    % rebind the loop (left leg is always left)
    sta = sta * LEFSystem.move_leg(loop_idx, new_pos1);
    sta = sta * LEFSystem.move_leg(loop_idx+LEFSystem.N, new_pos2);
    
    % regenerate events for the loop itself and for its new neighbours
    regenerate_all_loop_events(LEFSystem, evheap, loop_idx);
    regenerate_neighbours(LEFSystem, evheap, new_pos1);
    regenerate_neighbours(LEFSystem, evheap, new_pos2);
    
else
    disp(strcat('event_idx assumed a forbidden value :', num2str(event_idx)))
    sta = 0;
end
end