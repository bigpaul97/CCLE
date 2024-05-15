function l_length = GetLoopLength2(loop_parent,Loop_order,l_foot,r_foot,l_point,r_point,N)
% l_point is the left anchor of loop_parent
% This function calculate the total length of the left-of-PT child loops of the current parent loop

foot_bb = [];
for j = 1:N
    if Loop_order(j,2) == Loop_order(loop_parent,2)+1 &&...
            l_foot(Loop_order(j,1)) >= l_point &&...%l_foot(Loop_order(loop_parent,1)) &&...
            r_foot(Loop_order(j,1)) <= r_point%Anc(1) % Find the possible nearest child loop of the loop_parent
        foot_bb = [foot_bb,Loop_order(j,1)];
    end
end
l_length = 0;
for k = 1:length(foot_bb)
    l_length = l_length+r_foot(foot_bb(k))-l_foot(foot_bb(k)); % the length of the child loops
end

end