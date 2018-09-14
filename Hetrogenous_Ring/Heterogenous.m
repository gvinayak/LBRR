clc;
clear all;
warning off;

% Network Entire Field dimension (suppose 600 X 600)
xm=600; 
ym=600;
exp = 1;    %Determines whether it is to expand the ring or not
con =0;     %Determines whether it is to contract the ring or not
% Grid Size (suppose 50 X 50)
grid_xm=50;
grid_ym=50;

n=350;                    % Total number of Sensor Nodes in the field = 300
sink_move=10;             % Sink Movement in every 20th round
block_xm=xm/grid_xm;      % Calculating total number of Grids on x-axis side (i.e. 12)
block_ym=ym/grid_ym;      % Calculating total number of Grids on y-axis side (i.e. 12)
t1=block_xm;
t2=block_ym;

comm_radius_type1=50;     % Communication range of Type-1 sensor nodes is 50 meter
comm_radius_type2=70;     % Communication range of Type-2 sensor nodes is 70 meter
comm_radius_type3=90;     % Communication range of Type-3 sensor nodes is 90 meter

num_grid=block_xm * block_ym;               % Calculating Total Number of Grids over entire Network Area (i.e. 144)
node_grid=floor(n/num_grid);                % Average number of nodes in each grid (i.e. 2)
node_grid_remain=n-(node_grid * num_grid);  % Remaining nodes after deployment of avg no. of nodes in each grid (i.e. 12)

ring_radius=60;          % Default ring radius is 150 meter
ring_radius1=150;
ring_radius2=220;

% Transmission and Reception Energy (Eelec=Etx=Erx)
ETX=50*0.000000001;
ERX=50*0.000000001;

%Transmit Amplifier Energy
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;

do=sqrt(Efs/Emp);  % Computation of do

PS=4000;       % Packet Size (in bits)
CS = 64;        % Control Packet Size (in bits)
EDA=5*0.000000001; % Data Aggregation Energy

E1=1;                    % Initial Energy of Sensor Node in Joule
Power = 30;
comm_radius=80;           % Communication range of nodes is 80 meter
rmax=21;                % Maximum number of rounds
a1=1;                     % alpha
b1=2;                     % beta
Packet_Loss = 0;          % 
%s = 10;                    % For the sink movement.
mo=0.30;                  % Percentage of nodes that are Type-2 (30%)
m1=0.60;                  % Percentage of nodes that are Type-3 (60%)
type2=n*mo;               % Total number of Type-2 nodes in the network (n*30%)
type3=type2*m1;           % Total number of Type-3 nodes in the network (type2*60%)
type1=n-(type2+type3);    % Total number of Type-1 nodes in the network

scale=4;

%%% Position Allotment %%%
for i=1:(n-node_grid_remain)
    S(i).xd=rand * grid_xm;
    S(i).yd=rand * grid_ym;
end 

i=1;
a=-(t1/2);
b=-((t2/2)-1);
d=-(t1/2);
e=-(t2/2);
for m1=1:t2
    for k=1:1:block_xm
        if a*grid_xm < xm/2
           for j=1:node_grid
               S(i).xd=S(i).xd + a*grid_xm;
               S(i).yd=S(i).yd + d*grid_ym;
               i=i+1;
           end
           a=a+1;
        end 
    end

    for k=1:1:block_ym-1
        if b*grid_ym < ym/2
           for j=1:node_grid
               S(i).xd=S(i).xd + e*grid_xm; 
               S(i).yd=S(i).yd + b*grid_ym; 
               i=i+1;
           end
           b=b+1;
        end 
    end
    a=m1-(t1/2);
    b=m1-((t2/2)-1);
    d=d+1;
    e=e+1;
    block_xm=block_xm-1;
    block_ym=block_ym-1;
end

for i=(n-node_grid_remain)+1:n
    S(i).xd=unifrnd(-xm/2, xm/2);
    S(i).yd=unifrnd(-ym/2, ym/2);
end

%  for i=1:n
%      plot(S(i).xd,S(i).yd,'o')
%      hold on
%  end
% % S(1).xd
S(n+1).xd=0;  % x-axis position of Network Center
S(n+1).yd=0;  % y-axis position of Network Center
%%% Position Allotment Completed %%%
for i = 1:n
    S(i).distance = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
end
%%% Energy allotment %%%
Suffule_ID=randperm(n,n);
for i=1:type2
       S(Suffule_ID(i)).E=E1*(1+a1);
       S(Suffule_ID(i)).Energy=1;
end
for i=type2+1:type2+type3
       S(Suffule_ID(i)).E=E1*(1+b1);
       S(Suffule_ID(i)).Energy=2;
       S(Suffule_ID(i)).Power = comm_radius_type3;
end
for i=type2+type3+1:n
       S(Suffule_ID(i)).E=E1;
       S(Suffule_ID(i)).Energy=0;
       S(Suffule_ID(i)).Power = comm_radius_type1;
end
%%% Energy allotment complete %%%

% Entire Energy of the system
Ent_energy=(E1*type1)+ ((E1*(1+a1))*type2) + ((E1*(1+b1))*type2);

flag_first_dead=0;      % Dead nodes are assumed to be zero at the beginning
flag_all_dead=0;
s=0;
%plot(S(n+1).xd,S(n+1).yd,'rX','MarkerSize', scale*2.8,'LineWidth',2); % Plotting Network Center
%hold on

% for i = 1:n
%     if (S(i).Energy == 2)
%         plot(S(i).xd,S(i).yd,'gp','MarkerSize', scale*2.5,'LineWidth',1.3); 
%         hold on;
%     elseif (S(i).Energy == 1)
%         plot(S(i).xd,S(i).yd,'b+','LineWidth',1.3);
%         hold on;
%     elseif (S(i).Energy == 0)
%         plot(S(i).xd,S(i).yd,'ko','LineWidth',1.3);
%         hold on;
%     end
% end

sink.xd=randi([-xm/2 xm/2]);   % x-cordinates of the Sink
sink.yd=randi([-ym/2 ym/2]);   % y-cordinates of the Sink
S(n+2).xd=sink.xd;             % Sink is a (n+2)th node, x-axis position of sink node
S(n+2).yd=sink.yd;             % Sink is a (n+2)th node, y-axis position of sink node

dead2=0;     % Number of dead nodes
dead_a2=0;   % Number of dead Type-2 Nodes
dead_s2=0;   % Number of dead Type-3 Nodes
dead_n2=0;   % Number of dead Type-1 Nodes

save('cordinate_het.mat','S');

% Beginning of the function
for r=1:rmax
    r
    if (r==s*sink_move)
        % Random Sink Location   
        sink.xd=randi([-xm/2 xm/2]);   % x-cordinates of the Sink
        sink.yd=randi([-ym/2 ym/2]);   % y-cordinates of the Sink
        S(n+2).xd=sink.xd;             % Sink is a (n+2)th node, x-axis position of sink node
        S(n+2).yd=sink.yd;             % Sink is a (n+2)th node, y-axis position of sink node
        s=s+1;
    end
    
    if(r==1)
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
            if(S(ring_cand(i)).Energy==2)
                    f_ring(size(f_ring,2)+1) = ring_cand(i);  
                    break
            end
        end

        % For making Type-2 as a first starting ring node if Type-3 not present as a ring candidate
        for i=1:size(ring_cand,2)
            if (size(f_ring,2)==0)
                if((S(ring_cand(i)).Energy)==1)
                        f_ring(size(f_ring,2)+1) = ring_cand(i); 
                        break
                end
            end
        end
        
        % For making Type-1 as a first starting ring node if Type-3 and type-2 both not present as a ring candidate
        for i=1:size(ring_cand,2)
            if (size(f_ring,2)==0)
                if((S(ring_cand(i)).Energy)==0)
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
        Tran = path_find(S, C(1), C(size(C,2)));
        for k = 1:size(Tran,2)
            C(size(C,2) + 1) = Tran(k);
        end
        C = unique(C, 'stable');
        C(size(C,2) + 1) = C(1);
        
%         for i = 1:(size(C,2) - 1)
%             x = [S(C(i)).xd,S(C(i+1)).xd];
%             y = [S(C(i)).yd,S(C(i+1)).yd];
%             plot(x,y);
%             hold on
%         end
        continue;
    end
        if (S(i).E<=0)
            dead2=dead2+1;
            S(i).E=0;
            if(S(i).ENERGY==0)
               dead_n2=dead_n2+1;
            end
            if(S(i).ENERGY==1)
               dead_a2=dead_a2+1;
            end
            if(S(i).ENERGY==2)
               dead_s2=dead_s2+1;
            end
           % hold on;    
        end
        if (ring_radius2>=300)
            con = 1;
            exp = 0;
        end
        if (exp == 1)&&(r ~= 0)
            f1 = @Het_Expand_Ring;
            C = f1(S, C);
            f1 = @Update_ring;
            S = f1(S, C);
            ring_radius1 = ring_radius1 + 40;
            ring_radius2 = ring_radius2 + 40;
        end
        if (ring_radius1<=40)
            con = 0;
            exp = 1;
        end
        if (con == 1)&&(r ~= 0)
            f1 = @Het_Contract_Ring;
            C = f1(S, C);
            f1 = @Update_ring;
            S = f1(S, C);
            ring_radius1 = ring_radius1 - 40;
            ring_radius2 = ring_radius2 - 40;
        end
        RN=size(C,2)-1; % Total number of Ring Nodes
        Anch = Det_Anch_Het(S, n, RN, comm_radius_type1, comm_radius_type2, comm_radius_type3);
        
%         for i = 1:(size(C,2) - 1)
%             x = [S(C(i)).xd,S(C(i+1)).xd];
%             y = [S(C(i)).yd,S(C(i+1)).yd];
%             plot(x,y,'r');
%             hold on
%         end
%       plot(S(n+2).xd,S(n+2).yd,'r^','MarkerSize', scale*2.5,'MarkerFaceColor',[1,0,0]); % Plotting Sink
%      for i = 1:size(Anch,2)
%         plot(S(Anch(i)).xd,S(Anch(i)).yd,'r*');
%         hold on;
%      end
     
     %%% Finding nearest ring node for each element
     f1 = @Find_ring;
     [S, c1] = f1(S, C);
     
     %%% Finding the nearest ring node for the first anchor
     f1 = @Find_Ring_Anch; 
     rng_ele = f1(S, C, Anch);
     
     %%% Finding the path
     f1 = @path_find;
     Rng_Anch_Path = f1(S, Anch(1), rng_ele);
     
     %%% Reducing energy due to control packet
        for i = 1:size(Rng_Anch_Path,2)-1
            distance = sqrt((S(Rng_Anch_Path(i)).xd - S(Rng_Anch_Path(i+1)).xd)^2 + (S(Rng_Anch_Path(i)).yd - S(Rng_Anch_Path(i+1)).yd)^2); 
            if (i == 1)
                if (distance >=do)
                    S(Rng_Anch_Path(i)).E = S(Rng_Anch_Path(i)).E - ((ETX+EDA)*(CS) + Emp * CS *(distance * distance * distance * distance)); 
                else
                    S(Rng_Anch_Path(i)).E = S(Rng_Anch_Path(i)).E - ((ETX+EDA)*(CS) + Emp * CS *(distance * distance)); 
                end
            end
            if (distance >=do)
                S(Rng_Anch_Path(i)).E = S(Rng_Anch_Path(i)).E - ((ETX+EDA)*(CS) + Emp * CS *(distance * distance * distance * distance)); 
                S(Rng_Anch_Path(i)).E = S(Rng_Anch_Path(i)).E - ((ERX+EDA)*(CS) + Emp * CS *(distance * distance * distance * distance)); 
            else
                S(Rng_Anch_Path(i)).E = S(Rng_Anch_Path(i)).E - ((ETX+EDA)*(CS) + Emp * CS *(distance * distance)); 
                S(Rng_Anch_Path(i)).E = S(Rng_Anch_Path(i)).E - ((ERX+EDA)*(CS) + Emp * CS *(distance * distance)); 
            end
        end
        if (distance >=do)
            S(Rng_Anch_Path(i + 1)).E = S(Rng_Anch_Path(i + 1)).E - ((ETX+EDA)*(CS) + Emp * CS *(distance * distance * distance * distance)); 
        else
            S(Rng_Anch_Path(i + 1)).E = S(Rng_Anch_Path(i + 1)).E - ((ETX+EDA)*(CS) + Emp * CS *(distance * distance)); 
        end
        f1 = @Find_Anch; 
        [S, c2] = f1(S, Anch);
        hop = [];
        hop(r) = c1 + c2;
        hop(r)
        %%% Transmission for each node towards its anchor node
        for i = 1:n
            f1 = @path_find;
            Anch_Path = f1(S, i, S(i).r);
            for j = 1:size(Anch_Path,2)-1
                distance = sqrt((S(Anch_Path(j)).xd - S(Anch_Path(j+1)).xd)^2 + (S(Anch_Path(j)).yd - S(Anch_Path(j+1)).yd)^2);
                if (i == 1)
                    if (distance >=do)
                        S(Anch_Path(j)).E = S(Anch_Path(j)).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
                    else
                        S(Anch_Path(j)).E = S(Anch_Path(j)).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance)); 
                    end
                end
                if (distance >=do)
                    S(Anch_Path(j)).E = S(Anch_Path(j)).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
                    S(Anch_Path(j)).E = S(Anch_Path(j)).E - ((ERX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
                else
                    S(Anch_Path(j)).E = S(Anch_Path(j)).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance)); 
                    S(Anch_Path(j)).E = S(Anch_Path(j)).E - ((ERX+EDA)*(PS) + Emp * PS *(distance * distance)); 
                end
            end
            if (distance >=do)
                S(Anch_Path(j + 1)).E = S(Anch_Path(j + 1)).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
            else
                S(Anch_Path(j + 1)).E = S(Anch_Path(j + 1)).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance)); 
            end
        end
        for i = 1: size(Anch, 2)
            if (S(Anch(i)).E <= 0.5)
                Packet_Loss = Packet_Loss + 1;
                Anch = Det_Anch_Het(S, n, RN, comm_radius_type1, comm_radius_type2, comm_radius_type3);
            end
        end
        
        NBS=50;                   % Node buffer size (100 packet)
        PSR= 150;                 % Packet service rate (150 packets/sec)
        AN=0;
        PAR=n-(RN+AN);
        K=PAR-PSR;
        AN=ceil(K/NBS);
        PAR=n-(RN+AN);
        Interference_Loss = ceil(0.2*PAR);
        Pack_Sent(1,r) = PAR;
        Pack_Loss(1,r) = Packet_Loss+Interference_Loss;

        RE=([S.E].');
        OV_RE(1,r)=sum(RE);
        Dead(1,r)=dead2;
        
end 
Sent=sum(Pack_Sent);
Loss=sum(Pack_Loss);
Packet_drop_ratio5= Loss/Sent;

Residual5=sum(OV_RE);
Energy_Consumption5=Ent_energy-Residual5;

OV_dead5=sum(Dead);
OV_Alive5=n-OV_dead5;