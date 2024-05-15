function loop_order = LoopinLoop(l_foot,r_foot,N)

% determine whether a loop is a parent loop
% loop_order: 1st column, the loop index
% loop_order: 2nd column, indicates the number of parent loops
% loop_order: 3rd column, the index of the closest parent loop

loop_order = zeros(N,3);

for i = 1:N
    parents = 0;
    parents_lfoot = 0;
    
    loop_order(i,1) = i;
    for j = 1:N
        if l_foot(i) > l_foot(j) && r_foot(i) < r_foot(j)
            loop_order(i,2) = loop_order(i,2)+1;
            if l_foot(j) > parents_lfoot
                parents_lfoot = l_foot(j);
                parents = j;
            end
        end
    end
    
    loop_order(i,3) = parents;
    
end

end
