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
sink_move=20;             % Sink Movement in every 20th round
block_xm=xm/grid_xm;      % Calculating total number of Grids on x-axis side (i.e. 12)
block_ym=ym/grid_ym;      % Calculating total number of Grids on y-axis side (i.e. 12)
t1=block_xm;
t2=block_ym;
dead = 0;

% Transmission and Reception Energy (Eelec=Etx=Erx)
ETX=50*0.000000001;
ERX=50*0.000000001;

%Transmit Amplifier Energy
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
do=sqrt(Efs/Emp);  % Computation of do

PS=4000;       % Packet Size (in bits)
CS = 32;        % Control Packet Size (in bits)
EDA=5*0.000000001; % Data Aggregation Energy

num_grid=block_xm * block_ym;               % Calculating Total Number of Grids over entire Network Area (i.e. 144)
node_grid=floor(n/num_grid);                % Average number of nodes in each grid (i.e. 2)
node_grid_remain=n-(node_grid * num_grid);  % Remaining nodes after deployment of avg no. of nodes in each grid (i.e. 12)

ring_radius=60;          % Default ring radius is 150 meter
ring_radius1=150;
ring_radius2=220;
Eo=10;                    % Initial Energy of Sensor Node in Joule
comm_radius=80;           % Communication range of nodes is 80 meter
rmax=20;                  % Maximum number of rounds
s = 1;                    % For the sink movement.
for i=1:(n-node_grid_remain)
    S(i).xd=rand * grid_xm;
    S(i).yd=rand * grid_ym;
end 

i=1;
a=-(t1/2);
b=-((t2/2)-1);
d=-(t1/2);
e=-(t2/2);
for m=1:t2
    
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
    a=m-(t1/2);
    b=m-((t2/2)-1);
    d=d+1;
    e=e+1;
    block_xm=block_xm-1;
    block_ym=block_ym-1;
end

for i=(n-node_grid_remain)+1:n
    S(i).xd=unifrnd(-xm/2, xm/2);
    S(i).yd=unifrnd(-ym/2, ym/2);
end


% S(1).xd
S(n+1).xd=0;  % x-axis position of Network Center
S(n+1).yd=0;  % y-axis position of Network Center
flag_first_dead=0;      % Dead nodes are assumed to be zero at the beginning
flag_all_dead=0;
s=0;
Packet_Loss = 0;
% plot(S(n+1).xd,S(n+1).yd,'kX'); % Plotting Network Center
% hold on
for i = 1:n
    S(i).type = 'reg';   % Initially all nodes are Regular Nodes
    S(i).Energy=Eo;     
end
sink.xd=randi([-xm/2 xm/2]);   % x-cordinates of the Sink
sink.yd=randi([-ym/2 ym/2]);   % y-cordinates of the Sink
S(n+2).xd=sink.xd;             % Sink is a (n+2)th node, x-axis position of sink node
S(n+2).yd=sink.yd;             % Sink is a (n+2)th node, y-axis position of sink node

% Beginning of the function
for r=1:rmax
     if (r==s*sink_move)
        % Random Sink Location   
        sink.xd=randi([-xm/2 xm/2]);   % x-cordinates of the Sink
        sink.yd=randi([-ym/2 ym/2]);   % y-cordinates of the Sink
        S(n+2).xd=sink.xd;             % Sink is a (n+2)th node, x-axis position of sink node
        S(n+2).yd=sink.yd;             % Sink is a (n+2)th node, y-axis position of sink node
        s=s+1;
     end
%      plot(S(n+2).xd,S(n+2).yd,'g^'); % Plotting Sink
     hold on;
     if (r == 1)
        f1 = @Iden_ring;
        S = f1(S, n, ring_radius1, ring_radius2);     % Potential Ring Elements
        f1 = @Ourring;
        C = f1(S);                                  % Determine Ring Elements
        C = Ourring(S);
        continue;
     end
     if (S(i).Energy<=0)
            dead = dead+1;
            S(i).E=0;
           % hold on;    
     end
     if (ring_radius2>=300)
        con = 1;
        exp = 0;
     end
     if (exp == 1)&&(r ~= 1)
        f1 = @Expand_Ring;
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
        f1 = @Contract_Ring;
        C = f1(S, C);
        f1 = @Update_ring;
        S = f1(S, C);
        ring_radius1 = ring_radius1 - 40;
        ring_radius2 = ring_radius2 - 40;
     end
     
%      Sink selects the nearest node as anchor node
     f1 = @Det_Anch;
     anch = f1(S);
     
     for i = 1:n
     % Source gets the anchor node's position from the ring
     % Determining the closest ring node for each node
        f1 = @Find_ring;
        [S, c1] = f1(S, C);
        
     % For that it needs the path from source to the nearest ring element
        f1 = @path_find;
        Tran_r = f1(S, i, S(i).r);
        
        for j = 1:size(Tran_r,2)-1
            distance = sqrt((S(Tran_r(j)).xd - S(Tran_r(j+1)).xd)^2 + (S(Tran_r(j)).yd - S(Tran_r(j+1)).yd)^2); 
            if (j == 1)
                if (distance >=do)
                    S(Tran_r(j)).Energy = S(Tran_r(j)).Energy - ((ETX+EDA)*(CS) + Efs * CS *(distance * distance * distance * distance)); 
                else
                    S(Tran_r(j)).Energy = S(Tran_r(j)).Energy - ((ETX+EDA)*(CS) + Efs * CS *(distance * distance)); 
                end
            end
            if (distance >=do)
                S(Tran_r(j)).Energy = S(Tran_r(j)).Energy - ((ETX+EDA)*(CS) + Efs * CS *(distance * distance * distance * distance)); 
                S(Tran_r(j)).Energy = S(Tran_r(j)).Energy - ((ERX+EDA)*(CS) + Efs * CS *(distance * distance * distance * distance)); 
            else
                S(Tran_r(j)).Energy = S(Tran_r(j)).Energy - ((ETX+EDA)*(CS) + Efs * CS *(distance * distance)); 
                S(Tran_r(j)).Energy = S(Tran_r(j)).Energy - ((ERX+EDA)*(CS) + Efs * CS *(distance * distance)); 
            end
        end
        if (distance >=do)
            S(Tran_r(j + 1)).Energy = S(Tran_r(j + 1)).Energy - ((ETX+EDA)*(CS) + Efs * CS *(distance * distance * distance * distance)); 
        else
            S(Tran_r(j + 1)).Energy = S(Tran_r(j + 1)).Energy - ((ETX+EDA)*(CS) + Efs * CS *(distance * distance)); 
        end
        
        % Source transmits to the anchor node
        c2 = [];
        f1 = @path_find;
        Tran_a = f1(S, i, anch);
        c2(i) = size(Tran_a, 2); 
        
        for j = 1:size(Tran_a,2)-1
            distance = sqrt((S(Tran_a(j)).xd - S(Tran_a(j+1)).xd)^2 + (S(Tran_a(j)).yd - S(Tran_a(j+1)).yd)^2); 
            if (j == 1)
                if (distance >=do)
                    S(Tran_a(j)).Energy = S(Tran_a(j)).Energy - ((ETX+EDA)*(PS) + Efs * PS *(distance * distance * distance * distance)); 
                else
                    S(Tran_a(j)).Energy = S(Tran_a(j)).Energy - ((ETX+EDA)*(PS) + Efs * PS *(distance * distance)); 
                end
            end
            if (distance >=do)
                S(Tran_a(j)).Energy = S(Tran_a(j)).Energy - ((ETX+EDA)*(PS) + Efs * PS *(distance * distance * distance * distance)); 
                S(Tran_a(j)).Energy = S(Tran_a(j)).Energy - ((ERX+EDA)*(PS) + Efs * PS *(distance * distance * distance * distance)); 
            else
                S(Tran_a(j)).Energy = S(Tran_a(j)).Energy - ((ETX+EDA)*(PS) + Efs * PS *(distance * distance)); 
                S(Tran_a(j)).Energy = S(Tran_a(j)).Energy - ((ERX+EDA)*(PS) + Efs * PS *(distance * distance)); 
            end
        end
        if (distance >=do)
            S(Tran_a(j + 1)).Energy = S(Tran_a(j + 1)).Energy - ((ETX+EDA)*(PS) + Efs * PS *(distance * distance * distance * distance)); 
        else
            S(Tran_a(j + 1)).Energy = S(Tran_a(j + 1)).Energy - ((ETX+EDA)*(PS) + Efs * PS *(distance * distance)); 
        end
     end
     %%%%%%%%%%%%%%%%%%%%%%%EDIT here
     c22 = sum(c2);
     hop(r) = c1 + c22;
     hop(r)
     %S(anch).Energy
     if (S(anch).Energy <= 1)
        Packet_Loss = Packet_Loss + 1;
        anch = Det_Anch(S);
     end
     S(r).Packet_Loss = Packet_Loss;
end
% for i = 1:rmax
%     S(r).Packet_Loss
% end