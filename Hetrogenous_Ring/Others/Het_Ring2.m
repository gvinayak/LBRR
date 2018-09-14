function [ C ] = Het_Ring2( S )
% This function determines nodes for the ring.
    %For selecting first ring node which has the minimum distance from Network Center 
    temp=9999; 
    n = 300;
    R = [];
    C = [];
    for i=1:n
        if(strcmp(S(i).type,'pot_ring')==1)
            if (S(i).distance < temp)
                temp = S(i).distance;
                f_ring = i;
            end
        end
    end
    plot(S(f_ring).xd,S(f_ring).yd,'yellow .'); % Plotting first ring node
    f = 0;
    for i = 1:n
        if(strcmp(S(i).type,'pot_ring')==1)
            f = f + 1;
            R(f).id = i; 
            if(S(i).ang >= S(f_ring).ang)
                R(f).ang = S(i).ang - S(f_ring).ang;
            else
                R(f).ang = 360 - S(f_ring).ang + S(i).ang;
            end
        end
    end
    dist = 80;
    C(1) = f_ring;
    current = C(1);
    r_end = 0;
    tr = 0;
% % %         Ring_Construction
        while (r_end~=1)
            Tran = [];
            for i = 1:size(R,2)
                if(current == R(i).id)
                  cur_id = i;
                  break
                end
            end
            max_ang = 0;
            f = 0;
            for i=1:size(R,2)
                d = sqrt((S(R(i).id).xd - S(current).xd)^2 + (S(R(i).id).yd - S(current).yd)^2);
                wid = R(i).ang - R(cur_id).ang;
                if(d <= dist)
                    if (R(i).id == C(1))&&(wid< -300)
                        r_end = 1;
                        break;
                    elseif(wid<=90)&&(wid>=max_ang)
                        max_ang = wid;
                        maxid = R(i).id;
                        f = 1;
                    end
                end
            end
            if(tr >= 360)||(r_end == 1)
                r_end = 1;
                break
            end
            tr = tr + max_ang;
            current = maxid;
            C(size(C,2) + 1) = current;
%             r_end = 1;
        end
        Tran = path(S, C(1), C(size(C,2)));
        for k = 1:size(Tran,2)
            C(size(C,2) + 1) = Tran(k);
        end
        C = unique(C, 'stable');
        C(size(C,2) + 1) = C(1);
%     for i = 1:size(C,2)
%         plot(S(C(i)).xd,S(C(i)).yd,'*-');
%     end
    for i = 1:(size(C,2) - 1)
        x = [S(C(i)).xd,S(C(i+1)).xd];
        y = [S(C(i)).yd,S(C(i+1)).yd];
        plot(x,y);
    end
%     C = Expand_Ring(S, C);
end    
