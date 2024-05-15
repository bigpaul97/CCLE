function Cor_NP = SpatialCor_TwoPoints(l_foot,r_foot,PT1,PT2,N)

if PT1>PT2
    PT3=PT1;
    PT1=PT2;
    PT2=PT3;
end

Loop_order = LoopinLoop(l_foot,r_foot,N);

parent1 = GetParentLoop(l_foot,r_foot,PT1,N); % nearest parent of PT1
parent2 = GetParentLoop(l_foot,r_foot,PT2,N); % nearest parent of PT2

NP = [];
if parent1 == 0 && parent2 == 0
    NP = [NP,GetBackboneLength2(l_foot,r_foot,PT1,PT2,N)];
else 
    if (parent1 == 0) + (parent2 == 0) == 1
        if parent1 == 0
            [NP,Anc] = Calculate_NP(parent2,Loop_order,l_foot,r_foot,PT2,NP,N);
            NP = [NP,GetBackboneLength2(l_foot,r_foot,PT1,Anc(1),N)];
        else
            [NP,Anc] = Calculate_NP(parent1,Loop_order,l_foot,r_foot,PT1,NP,N);
            NP = [NP,GetBackboneLength2(l_foot,r_foot,Anc(2),PT2,N)];
        end
    else
        loop_parent1 = parent1;
        loop_parent2 = parent2;
        Anc1(1) = PT1;
        Anc1(2) = PT1;
        Anc2(1) = PT2;
        Anc2(2) = PT2;
        while Loop_order(loop_parent1,2) ~= Loop_order(loop_parent2,2)
            if Loop_order(loop_parent1,2) > Loop_order(loop_parent2,2)
                [NP,Anc1,loop_parent1] = Calculate_SingleNP(loop_parent1,Loop_order,l_foot,r_foot,Anc1,NP,N);
            else
                [NP,Anc2,loop_parent2] = Calculate_SingleNP(loop_parent2,Loop_order,l_foot,r_foot,Anc2,NP,N);
            end
        end
        if Loop_order(loop_parent1,1) == Loop_order(loop_parent2,1)
            NP1 = Anc2(1)-Anc1(2)-GetLoopLength2(loop_parent1,Loop_order,l_foot,r_foot,Anc1(2),Anc2(1),N);
            NP2 = r_foot(loop_parent1)-l_foot(loop_parent1)-GetLoopLength2(loop_parent1,Loop_order,l_foot,r_foot,l_foot(loop_parent1),r_foot(loop_parent1),N)-NP1;
%             if NP1 == 0
%                 NP1 =1;
%             end
%              if NP2 == 0
%                 NP2 =1;
%             end
            NP = [NP,NP1,NP2];
        else
            while loop_parent1 ~= loop_parent2
                [NP,Anc1,loop_parent1] = Calculate_SingleNP(loop_parent1,Loop_order,l_foot,r_foot,Anc1,NP,N);
                [NP,Anc2,loop_parent2] = Calculate_SingleNP(loop_parent2,Loop_order,l_foot,r_foot,Anc2,NP,N);
            end
            if loop_parent1 ~= 0
                NP1 = Anc2(1)-Anc1(2)-GetLoopLength2(loop_parent1,Loop_order,l_foot,r_foot,Anc1(2),Anc2(1),N);
                NP2 = r_foot(loop_parent1)-l_foot(loop_parent1)-GetLoopLength2(loop_parent1,Loop_order,l_foot,r_foot,l_foot(loop_parent1),r_foot(loop_parent1),N)-NP1;
                NP = [NP,NP1,NP2];
            else
                NP = [NP,GetBackboneLength2(l_foot,r_foot,Anc1(2),Anc2(1),N)];
            end
        end
    end
end

Cor_NP = 0;
for i = 1:2:length(NP)-1
%     if NP(i)==0
%         NP(i)=1;
%     end
%     if NP(i+1)==0
%         NP(i+1)=1;
%     end
    Cor_NP = Cor_NP+NP(i)*NP(i+1)/(NP(i)+NP(i+1));
end

if mod(length(NP),2) == 1
    Cor_NP = Cor_NP+NP(end);
end
