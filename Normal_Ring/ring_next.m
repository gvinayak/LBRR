function [ next, r_end, max_ang] = ring_next(S, R, current, cur_id, f_ring)
%Finds the next element of the ring
    max_ang2 = 0;max_ang1 = 0;max_ang0 = 0;
    r_end = 0;
    found2 = 0;found1 = 0;found0 = 0;
    for i=1:size(R,2)
        d = sqrt((S(R(i).id).xd - S(current).xd)^2 + (S(R(i).id).yd - S(current).yd)^2);
        wid = R(i).ang - R(cur_id).ang;
        if(R(i).type == 2)
            if(d <= R(cur_id).Pow + R(i).Pow)
                if (R(i).id == f_ring)&&(wid > 0)
                    r_end = 1;
                    break;
                elseif(wid<=90)&&(wid>=max_ang2)
                    max_ang2 = wid;
                    maxid2 = R(i).id;
                    found2=1;
                end
            end
        elseif(R(i).type == 1)
            if(d <= R(cur_id).Pow + R(i).Pow)||(found2~=0)
                if (R(i).id == f_ring)&&(wid > 0)
                    r_end = 1;
                    break;
                elseif(wid<=90)&&(wid>=max_ang1)
                    max_ang1 = wid;
                    maxid1 = R(i).id;
                    found1=0;
                end
            end
        elseif(R(i).type == 0)
            if(d <= R(cur_id).Pow + R(i).Pow)||(found2~=0)||(found1~=0)
                if (R(i).id == f_ring)&&(wid > 0)
                    r_end = 1;
                    break;
                elseif(wid<=90)&&(wid>=max_ang0)
                    max_ang0 = wid;
                    maxid0 = R(i).id;
                    found0=1;
                end
            end
        end
    end
    if (r_end ==1)
        next = f_ring;
        max_ang = 0;
    elseif (found2==1)
        next = maxid2;
        max_ang = max_ang2;
    elseif(found2==0)&&(found1==1)
        next = maxid1;
        max_ang = max_ang1;
    elseif(found2==0)&&(found1==0)&&(found0==1)
        next = maxid0;
        max_ang = max_ang0;
    end
end