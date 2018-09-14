function [ Anch ] = Het_One_Anch( S )
%Find one anchor node
    Max_E = 0;
    n = 350;
    mind = 999;
    for i = 1:n
        d = sqrt((S(i).xd - S(n+2).xd)^2 + (S(i).yd - S(n+2).yd)^2);
        if (d<= S(i).Power)&&(mind > d) 
            if (Max_E <= S(i).E)
                Max_E = d;
                mind = d;
                Anch = i;
            end
        end
    end
end