function [ S ] = Iden_ring( S, n, ring_radius1, ring_radius2)
%This function determines the potential ring members.
for i = 1:n
    distance = sqrt(((S(i).xd - S(n+1).xd)^2) + ((S(i).yd - S(n+1).yd)^2));
    S(i).distance = distance;
    ang = (180/pi)*atan2(S(i).yd, S(i).xd);
    if ang<0
        ang = ang + 360;
    end
    S(i).ang = ang;
    if (distance>=ring_radius1 && distance<=ring_radius2)
        S(i).type='pot_ring';
    else
        S(i).type='reg';
    end
end
end