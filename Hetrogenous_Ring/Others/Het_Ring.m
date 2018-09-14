function [ C ] = Het_Ring( S )
% This function determines nodes for the ring.
    %For selecting first ring node which has the minimum distance from Network Center 
    n = 350;
%     if(r==1)
        % Determining initial ring candidates that are near the initial ring radius
        ring_cand=[];
        for i = 1:n
            distance = sqrt(((S(i).xd - S(n+1).xd)^2) + ((S(i).yd - S(n+1).yd)^2));
            S(i).distance = distance;                % Every node's distance from network center (NC) is stored
            ang = (180/pi)*atan2(S(i).yd, S(i).xd);  % Calculating every node's angle w.r.t. NC
            if ang<0
               ang = ang + 360;         % If angle is negative then becoming it positve so that calculted angle will be from positive x-axis
            end
            S(i).ang = ang;             % Every node's angle from positive x-axis is stored in variable S(i).ang
            if (distance>=ring_radius1 && distance<=ring_radius2)
                S(i).role='ring_cand';
                ring_cand(size(ring_cand,2)+1)=i;
            end
        end

        % For making Type-3 as a first starting ring node
        f_ring=[];
        for i=1:size(ring_cand,2)
            if(strcmp(S(ring_cand(i)).Energy)==2)
                    f_ring(size(f_ring,2)+1) = ring_cand(i);  
                    break
            end
        end

        % For making Type-2 as a first starting ring node if Type-3 not present as a ring candidate
        for i=1:size(ring_cand,2)
            if (size(f_ring,2)==0)
                if(strcmp(S(ring_cand(i)).Energy)==1)
                        f_ring(size(f_ring,2)+1) = ring_cand(i); 
                        break
                end
            end
        end
        
        % For making Type-1 as a first starting ring node if Type-3 and type-2 both not present as a ring candidate
        for i=1:size(ring_cand,2)
            if (size(f_ring,2)==0)
                if(strcmp(S(ring_cand(i)).Energy)==0)
                        f_ring(size(f_ring,2)+1) = ring_cand(i);  
                        break
                end
            end
        end
        
        % For calculating angle between first ring node and other ring candidates 
        C=[];
       
        for i = 1:size(ring_cand,2)
            if(strcmp(S(ring_cand(i)).role,'ring_cand')==1)
               if(S(ring_cand(i)).ang >= S(f_ring(1)).ang)
                  R(i).ang = S(ring_cand(i)).ang - S(f_ring(1)).ang;       % Every Potential Ring node's angle w.r.t. first ring node's angle (Difference between first rind node's angle and potential ring node's angle)
               else
                  R(i).ang = 360 - S(f_ring(1)).ang + S(ring_cand(i)).ang; % Every Potential Ring node's angle w.r.t. first ring node's angle
               end
            end
        end
        
        C(1) = f_ring(1);
        current = C(1);
        r_end = 0;
        tr = 0;
    
        % For Ring_Construction
        while (r_end~=1)
               Tran = [];
               for i = 1:size(ring_cand,2)
                   if(current == ring_cand(i))
                      cur_id = i;
                      break
                   end
               end
               max_ang = 0;
               f = 0;
               
               for i=1:size(ring_cand,2)
                   d = sqrt((S(ring_cand(i)).xd - S(current).xd)^2 + (S(ring_cand(i)).yd - S(current).yd)^2); 
                   wid = R(i).ang - R(cur_id).ang;
                   if(d < 2*comm_radius_type2)
                      if (ring_cand(i) == C(1))&&(wid< -300)
                          r_end = 1;
                          break;
                      elseif(wid<=90)&&(wid>=max_ang)
                             max_ang = wid;
                             maxid = ring_cand(i);
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
        end
        Tran = path(S, C(1), C(size(C,2)));
        for k = 1:size(Tran,2)
            C(size(C,2) + 1) = Tran(k);
        end
        C = unique(C, 'stable');
        C(size(C,2) + 1) = C(1);
        for i = 1:size(C,2)
            %plot(S(C(i)).xd,S(C(i)).yd,'ro');
        end
        for i = 1:(size(C,2) - 1)
            x = [S(C(i)).xd,S(C(i+1)).xd];
            y = [S(C(i)).yd,S(C(i+1)).yd];
            plot(x,y,'r');
            hold on
        end
%     end 
end    
