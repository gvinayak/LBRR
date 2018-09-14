function [ anch ] = Det_Anch( S )
%This function determines the nearestnode to the sink and makes it the
%anchor node.
    mind = 9999;
    n = 300;
    for i = 1:n
        d = sqrt((S(i).xd - S(n+2).xd)^2 + (S(i).yd - S(n+2).yd)^2);
        if d<=mind
            mind = d;
            anch = i;
        end
    end
end