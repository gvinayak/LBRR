function [ S, c2] = Find_Anch( S, Anch)
%This function determines the nearest anchor node element for each node
    minid = 0;
    c2 = 0;
    count = [];
    for i = 1:size(S,2) - 2
        count1 = 0;
        mind = 9999;
        for j = 1:size(Anch,2)
            if(i == Anch(j))
                S(i).r = 0;
                continue
            else
                d = sqrt((S(Anch(j)).xd-(S(i).xd))^2 + (S(Anch(j)).yd-(S(i).yd))^2);
                if(d<mind)&&(d>0)
                    mind = d;
                    minid = Anch(j);
                    count1 = count1 + 1;
                end
            end
        end
        S(i).r = minid;
        count(i) = count1;
    end
    c2 = sum(count);
end