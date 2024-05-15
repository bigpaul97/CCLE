function bb_length = GetBackboneLength2(l_foot,r_foot,l_point,r_point,N)

loop_order = LoopinLoop(l_foot,r_foot,N);

foot_bb = [];
for j = 1:N
    if loop_order(j,2) == 0 && l_foot(loop_order(j,1)) >= l_point &&...
            r_foot(loop_order(j,1)) <= r_point
        foot_bb = [foot_bb,loop_order(j,1)];
    end
end

loop_length = 0;
for k = 1:length(foot_bb)
    loop_length = loop_length+r_foot(foot_bb(k))-l_foot(foot_bb(k));
end
bb_length = r_point-l_point-loop_length;

% if bb_length == 0
%     bb_length =1;
% end