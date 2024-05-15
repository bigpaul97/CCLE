function [NP,Anc] = Calculate_NP(loop_parent,Loop_order,l_foot,r_foot,PT,NP,N)

Anc(1) = PT;
Anc(2) = PT;
while loop_parent
    [NP,Anc,loop_parent] = Calculate_SingleNP(loop_parent,Loop_order,l_foot,r_foot,Anc,NP,N);
end

end