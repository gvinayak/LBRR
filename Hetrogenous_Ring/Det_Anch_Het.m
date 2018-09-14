% Sink selects multiple agent nodes
function [ AN_ID ] = Det_Anch_Het( S, n, RN, comm_radius_type1, comm_radius_type2, comm_radius_type3)
    NBS=50;                   % Node buffer size (100 packet)
    PSR= 150;                 % Packet service rate (150 packets/sec)
    AN=0;
    PAR=n-(RN+AN);
    K=PAR-PSR;
    AN=ceil(K/NBS);

    Sink_ANRange=[];
    for i=1:n
        sink_dist=sqrt(((S(i).xd - S(n+2).xd)^2) + ((S(i).yd - S(n+2).yd)^2));
        if (S(i).Energy==2)
            if (sink_dist<comm_radius_type3)
                Sink_ANRange(size(Sink_ANRange,2)+1)=i;
            end
        end
        if (S(i).Energy==1)
            if (sink_dist<comm_radius_type2)
                Sink_ANRange(size(Sink_ANRange,2)+1)=i;
            end
        end
        if (S(i).Energy==0)
            if (sink_dist<comm_radius_type1)
                Sink_ANRange(size(Sink_ANRange,2)+1)=i;
            end
        end
    end
    
    AN_ID=[];
    for i=1:size(Sink_ANRange,2)
        if (S(Sink_ANRange(i)).Energy==2)
            AN_ID(size(AN_ID,2)+1)=Sink_ANRange(i);
        end
        if (size(AN_ID,2)==AN)
           break
        end  
    end
    
    if (size(AN_ID,2)<AN)
        for i=1:size(Sink_ANRange,2)
            if (S(Sink_ANRange(i)).Energy==1)
                AN_ID(size(AN_ID,2)+1)=Sink_ANRange(i);
            end
            if (size(AN_ID,2)==AN)
               break
            end  
        end
    end
    
    if (size(AN_ID,2)<AN)
        for i=1:size(Sink_ANRange,2)
            if (S(Sink_ANRange(i)).Energy==0)
                AN_ID(size(AN_ID,2)+1)=Sink_ANRange(i);
            end
            if (size(AN_ID,2)==AN)
               break
            end  
        end
    end
    
%     for i = 1:size(AN_ID,2)
%         plot(S(AN_ID(i)).xd,S(AN_ID(i)).yd,'k*');
%         hold on
%     end
end