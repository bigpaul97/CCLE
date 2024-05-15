function regenerate_event(LEFSystem, evheap, event_idx)

%     Regenerate an event in an event heap. If the event is currently impossible (e.g. a step
%     onto an occupied site) then the new event is not created, but the existing event is not
%     modified.
%
%     Possible events:
%     1 to 2N : a step to the left
%     2N+1 to 4N : a step to the right
%     4N+1 to 5N : passive unbinding
%     5N+1 to 6N : rebinding to a randomly chosen site


% A step to the left or to the right.
if event_idx <= 4 * LEFSystem.N
    if event_idx <= 2 * LEFSystem.N
        leg_idx = event_idx;
        direction = -1;
    else
        leg_idx = event_idx - 2 * LEFSystem.N;
        direction = 1;
    end
    
    if LEFSystem.locs(leg_idx) > 0
        idxx1 = LEFSystem.locs(leg_idx);
        idxx2 = int64((3-direction)/2);
        % Local velocity = velocity * permeability
        local_vel = LEFSystem.perms(idxx1, idxx2)...
            * LEFSystem.vels(leg_idx + (direction+1) * LEFSystem.N);
%         disp(idxx1)
        if local_vel > 0
            new_position = LEFSystem.locs(leg_idx) + direction;
            if new_position == 0
                new_position = LEFSystem.L;
            elseif new_position > LEFSystem.L
                new_position = mod(new_position, LEFSystem.L);
            end

            if LEFSystem.lattice(new_position) <= 0
                evheap.add_event(event_idx,LEFSystem.time + exprnd(1/local_vel));
            end
        end
    end
    
    % Passive unbinding.
elseif event_idx > 4 * LEFSystem.N && event_idx <= 5 * LEFSystem.N
    loop_idx = event_idx - 4 * LEFSystem.N;
    if LEFSystem.locs(loop_idx) > 0 && LEFSystem.locs(loop_idx+LEFSystem.N) > 0
        evheap.add_event(event_idx,LEFSystem.time + exprnd(LEFSystem.lifespans(loop_idx)))
    end
    
    % Rebinding from the solution to a random site.
elseif event_idx > 5 * LEFSystem.N && event_idx <= 6 * LEFSystem.N
    loop_idx = event_idx - 5 * LEFSystem.N;
    if LEFSystem.locs(loop_idx) <= 0 && LEFSystem.locs(loop_idx+LEFSystem.N) <= 0
        evheap.add_event(event_idx,LEFSystem.time + exprnd(LEFSystem.rebinding_times(loop_idx)))
    end
end

end

