function [ S, c1 ] = Find_ring( S, C )
%This function determines the nearest ring node element for each node
    minid = 0;
    c1 = 0;
    count = [];
    for i = 1:size(S,2) - 2
        mind = 9999;
        count1 = 0;
        for j = 1:size(C,2)
            if(i == C(j))
                S(i).r = 0;
                continue
            else
                d = sqrt((S(C(j)).xd-(S(i).xd))^2 + (S(C(j)).yd-(S(i).yd))^2);
                if(d<mind)&&(d>0)
                    mind = d;
                    minid = C(j);
                    count1 = count1 + 1;
                end
            end
        end
        S(i).r = minid;
        count(i) = count1;
    end
    c1 = sum(count);
end