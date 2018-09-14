function [ f_ring, R ] = Het_Nearest( S )
%Find the nearest element for a heterogenous network, prioroty wise 
    f = 0;
    R = [];
    n = 350;
    found2 = 0;found1 = 0;found0 = 0;
    temp2 = 9999; temp1 = 9999; temp0 = 9999;
    for i = 1:n
        if(strcmp(S(i).type,'pot_ring')==1)
            f = f + 1;
            R(f).id = i;
            R(f).type = S(i).Energy;
            R(f).Pow = S(i).Power;
            R(f).ang = S(i).ang;
            if(S(i).Energy == 2)
                found2 = 1;
                if (temp2 > S(i).distance)
                    index2 = i;
                    temp2 = S(i).distance;
                end
            elseif(S(i).Energy == 1)
                found1 = 1;
                if(temp1 > S(i).distance)
                    index1 = i;
                    temp1 = S(i).distance;
                end
            elseif(S(i).Energy == 0)
                found0 = 1;
                if(temp0 > S(i).distance)
                    index0 = i;
                    temp0 = S(i).distance;
                end
            end
        end
    end
    if found2==1
        f_ring = index2;
    elseif(found2==0)&&(found1==1)
        f_ring = index1;
    elseif(found1==0)&&(found0==1)
        f_ring = index0;
    end
end