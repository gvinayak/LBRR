function [ rng_ele ] = Find_Ring_Anch(S, C, Anch)
%This function determines the nearest ring node element for each node
    minid = 0;
    mind = 9999;
        for j = 1:size(C,2)
            d = sqrt((S(C(j)).xd-(S(Anch(1)).xd))^2 + (S(C(j)).yd-(S(Anch(1)).yd))^2);
            if(d<mind)&&(d>0)
                mind = d;
                minid = C(j);
            end
        end
        rng_ele = minid;
end