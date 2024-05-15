function [NP,Anc,loop_parent] = Calculate_SingleNP(loop_parent,Loop_order,l_foot,r_foot,Anc,NP,N)

l_length = GetLoopLength2(loop_parent,Loop_order,l_foot,r_foot,l_foot(Loop_order(loop_parent,1)),Anc(1),N);
NP1 = Anc(1)-l_foot(Loop_order(loop_parent,1))-l_length;

r_length = GetLoopLength2(loop_parent,Loop_order,l_foot,r_foot,Anc(2),r_foot(Loop_order(loop_parent,1)),N);
NP2 = r_foot(Loop_order(loop_parent,1))-Anc(2)-r_length;

NP = [NP,NP1,NP2];
Anc(1) = l_foot(Loop_order(loop_parent,1));
Anc(2) = r_foot(Loop_order(loop_parent,1));
loop_parent = Loop_order(loop_parent,3);

end