

function [ ring ] = Det_Ring( S, C, anch )
%This function determines the nearestnode to the sink and makes it the
%anchor node.
    mind = 9999;
    for i = 1:size(C,2)
        d = sqrt((S(C(i)).xd - S(anch).xd)^2 + (S(C(i)).yd - S(anch).yd)^2);
        if d<mind
            mind = d;
            ring = i;
        end
    end
end