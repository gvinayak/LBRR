function [ next, ang ] = Comp_Ring( S, R, current )
%This function finds next ring element when not in range.
max_ang = 0;
maxd = 9999;
for i=1:size(S,2) - 2
    d = sqrt((S(i).xd - S(current).xd)^2 + (S(i).yd - S(current).yd)^2);
    wid = S(i).ang - R(current).ang;
    if(wid>max_ang)&&(d <= 80)
        minid = i;
        max_ang = wid;
    end
    if(strcmp(S(i).type,'pot_ring')==1)&&(wid>0)&&(d <= maxd)
        nt = i;
        maxd = d;
    end
end
next = path(S, minid, nt);
for i = 1:size(R,2)
    if(nt==R(i).id)
        ang = R(i).ang;
    end
end
    
end

