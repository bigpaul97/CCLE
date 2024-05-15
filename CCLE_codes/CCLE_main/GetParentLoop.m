function parent = GetParentLoop(l_foot,r_foot,PT,N)
% Return the nearest parent loop of the point PT

parent = 0;
parents_lfoot = 0;
for i = 1:N
    if l_foot(i) < PT && r_foot(i) > PT
        if l_foot(i) > parents_lfoot
            parents_lfoot = l_foot(i);
            parent = i;
        end
    end
end

