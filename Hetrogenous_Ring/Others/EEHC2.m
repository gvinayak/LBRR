clc;
clear all;
warning off;

% Field dimension
xm=100; 
ym=100; 

sink.x=xm/2;   % x-cordinates of the Sink
sink.y=ym/2;   % y-cordinates of the Sink
n=100;         % Number of Nodes in the field
PS=4000;       % Packet Size (in bits)
p=0.04;        % Optimal Election Probability of a node to become cluster head

% Energy Model (all values in Joules) 
E1=0.3;        % Initial Energy of node

% Transmission and Reception Energy (Eelec=Etx=Erx)
ETX=50*0.000000001;
ERX=50*0.000000001;

%Transmit Amplifier Energy
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;

EDA=5*0.000000001; % Data Aggregation Energy

%Values for Hetereogeneity
m=0.40;            % Percentage of nodes than are Type-2
u=0.60;            % Percentage of nodes than are Type-3
a=1;               % alpha
b=2;               % beta

rmax=9999;        % Maximum number of rounds
do=sqrt(Efs/Emp);  % Computation of do

%%%%%%%%%%%%%%%%%%%%%%%   EEHC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creation of the random Sensor Network
%figure(1);
for i=1:n
    S2(i).xd=randi([0 xm]);
    S2(i).yd=randi([0 ym]);
    S2(i).G=0;
    S2(i).type='N';  % Initially there are no cluster heads only Nodes
    
    % Energy Distribution for Type-1 Nodes
    if (i>= m*n + m*n*u + 1) 
        S2(i).E=E1;
        S2(i).ENERGY=0;
    end
    
    % Energy Distribution for Type-2 Nodes
    if (i< m*n + 1)  
        S2(i).E=E1*(1+a);
        S2(i).ENERGY=1;
    end
    
    % Energy Distribution for Type-3 Nodes
    if (i>= m*n + 1 && i< m*n + m*n*u + 1)
        S2(i).E=E1*(1+b);
        S2(i).ENERGY=2;
    end
end

S2(n+1).xd=sink.x;  % Sink is a (n+1)th node, x-axis position of sink node
S2(n+1).yd=sink.y;  % Sink is a (n+1)th node, y-axis position of sink node

rcountCHs2=0; % Counter for CHs per round
flag_first_dead2=0;
flag_all_dead2=0;

% Counter for bit transmitted to BS and to CHs
packets_TO_BS2=0;
packets_TO_CH2=0;

for r=0:rmax
    r
    % Election Probability for Type-1 Nodes
    pnrm=p/(1+(m*(a+(u*b))));
  
    % Election Probability for Type-2 Nodes
    padv=(p*(1+a))/(1+(m*(a+(u*b))));
    
    % Election Probability for Type-3 Nodes
    psup=(p*(1+b))/(1+(m*(a+(u*b))));
    
    % Operation for heterogeneous sub-epochs for Type-1 Nodes
    if(mod(r, round(1/pnrm))==0)
       for i=1:n
           if(S2(i).ENERGY==0)
              S2(i).G=0;
           end   
       end
    end

    % Operations for sub-epochs for Type-2 Nodes
    if(mod(r, round(1/padv))==0)
       for i=1:n
           if(S2(i).ENERGY==1)
              S2(i).G=0;
           end
       end
    end
    
    % Operation for heterogeneous sub-epochs for Type-3 Nodes
    if(mod(r, round(1/psup))==0)
       for i=1:n
           if(S2(i).ENERGY==2)
              S2(i).G=0;
           end   
       end
    end
    
    hold off;
    
    dead2=0;     % Number of dead nodes
    dead_a2=0;   % Number of dead Type-2 Nodes
    dead_s2=0;   % Number of dead Type-3 Nodes
    dead_n2=0;   % Number of dead Type-1 Nodes
    
    % Counter for bit transmitted to BS and to CHs per round
    PACKETS_TO_CH2(r+1)=0;
    PACKETS_TO_BS2(r+1)=0;

   % figure(1);

    for i=1:n
        %checking if there is a dead node
        if (S2(i).E<=0)
            %plot(S2(i).xd,S2(i).yd,'red .');
            dead2=dead2+1;
            S2(i).E=0;
            if(S2(i).ENERGY==0)
               dead_n2=dead_n2+1;
            end
            if(S2(i).ENERGY==1)
               dead_a2=dead_a2+1;
            end
            if(S2(i).ENERGY==2)
               dead_s2=dead_s2+1;
            end
           % hold on;    
        end
        if (S2(i).E>0)
            S2(i).type='N';  
            if (S2(i).ENERGY==0)  
              % plot(S2(i).xd,S2(i).yd,'bo');
            end
            if (S2(i).ENERGY==1)  
              % plot(S2(i).xd,S2(i).yd,'g+');
            end
            if (S2(i).ENERGY==2)  
             %  plot(S2(i).xd,S2(i).yd,'c<');
            end
          % hold on;
        end
    end
    
   % plot(S2(n+1).xd,S2(n+1).yd,'k^');

    DEAD2(r+1)=dead2;
    DEAD_N2(r+1)=dead_n2;
    DEAD_A2(r+1)=dead_a2;
    DEAD_S2(r+1)=dead_s2;
    ALIVE2(r+1)=n-dead2;

    % When the first node dies
    if (dead2==1)
        if(flag_first_dead2==0)
           first_dead2=r+1;
           flag_first_dead2=1;
        end
    end
    
    % When the first node dies
    if (dead2==n)
        if(flag_all_dead2==0)
           all_dead2=r+1;
           flag_all_dead2=1;
        end
    end

   countCHs2=0;  % Counter for CHs
   cluster2=1;   % No of clusters
   
   for i=1:n
       if(S2(i).E>0)
          temp_rand=rand;     
          if ((S2(i).G)==0)

              % Election of CHs in every round for Type-1 Nodes
              if((S2(i).ENERGY==0 && (temp_rand <= (pnrm/(1 - pnrm * mod(r,round(1/pnrm)))))))
                  countCHs2=countCHs2+1;
                  packets_TO_BS2=packets_TO_BS2+1;
                  S2(i).type='C';
                  S2(i).G=1;
                  C2(cluster2).xd=S2(i).xd;
                  C2(cluster2).yd=S2(i).yd;
                 % plot(S2(i).xd,S2(i).yd,'k*');
                  distance=sqrt((S2(i).xd-(S2(n+1).xd))^2 + (S2(i).yd-(S2(n+1).yd))^2);
                  C2(cluster2).distance=distance;
                  C2(cluster2).id=i;
                  cluster2=cluster2+1;
            
                  % Calculation of Energy dissipated  by CH nodes in case of Type-1 Nodes 
                  distance;
                  if (distance>do)
                      S2(i).E=S2(i).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
                  end
                  if (distance<=do)
                      S2(i).E=S2(i).E - ((ETX+EDA)*(PS)  + Efs * PS *(distance * distance)); 
                  end
              end     

              % Election of Cluster Heads for Type-2 nodes
              if((S2(i).ENERGY==1 && (temp_rand <= (padv/(1 - padv * mod(r,round(1/padv)))))))
                  countCHs2=countCHs2+1;
                  packets_TO_BS2=packets_TO_BS2+1;
                  S2(i).type='C';
                  S2(i).G=1;
                  C2(cluster2).xd=S2(i).xd;
                  C2(cluster2).yd=S2(i).yd;
                  %plot(S2(i).xd,S2(i).yd,'k*');
                  distance=sqrt((S2(i).xd-(S2(n+1).xd))^2 + (S2(i).yd-(S2(n+1).yd))^2);
                  C2(cluster2).distance=distance;
                  C2(cluster2).id=i;
                  cluster2=cluster2+1;
            
                 % Calculation of Energy dissipated  by CH nodes in case of Type-2 nodes
                 distance;
                 if (distance>do)
                     S2(i).E=S2(i).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
                 end
                 if (distance<=do)
                     S2(i).E=S2(i).E - ((ETX+EDA)*(PS)  + Efs * PS *(distance * distance)); 
                 end
              end
              
              % Election of CHs in every round for Type-3 Nodes
              if((S2(i).ENERGY==2 && (temp_rand <= (psup/(1 - psup * mod(r,round(1/psup)))))))
                  countCHs2=countCHs2+1;
                  packets_TO_BS2=packets_TO_BS2+1;
                  S2(i).type='C';
                  S2(i).G=1;
                  C2(cluster2).xd=S2(i).xd;
                  C2(cluster2).yd=S2(i).yd;
                %  plot(S2(i).xd,S2(i).yd,'k*');
                  distance=sqrt((S2(i).xd-(S2(n+1).xd))^2 + (S2(i).yd-(S2(n+1).yd))^2);
                  C2(cluster2).distance=distance;
                  C2(cluster2).id=i;
                  cluster2=cluster2+1;
            
                  % Calculation of Energy dissipated  by CH nodes in case of Type-1 Nodes 
                  distance;
                  if (distance>do)
                      S2(i).E=S2(i).E - ((ETX+EDA)*(PS) + Emp * PS *(distance * distance * distance * distance)); 
                  end
                  if (distance<=do)
                      S2(i).E=S2(i).E - ((ETX+EDA)*(PS)  + Efs * PS *(distance * distance)); 
                  end
              end   
              
          end
       end 
   end

   CLUSTERHEADS2(r+1)=cluster2-1;

   % Election of Associated  CH Nodes by each CM nodes
   for i=1:n
       if (S2(i).type=='N' && S2(i).E>0)
          if(cluster2-1>=1)
             min_dis=sqrt((S2(i).xd-S2(n+1).xd)^2 + (S2(i).yd-S2(n+1).yd)^2);
             min_dis_cluster=1;
             for c2=1:cluster2-1
                 temp2=min(min_dis, sqrt( (S2(i).xd-C2(c2).xd)^2 + (S2(i).yd-C2(c2).yd)^2));
                 if ( temp2 < min_dis )
                     min_dis=temp2;
                     min_dis_cluster=c2;
                 end
             end
       
            % Energy dissipated by associated CM nodes
            min_dis;
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(PS) + Emp * PS *( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(PS) + Efs * PS *( min_dis * min_dis)); 
            end
            
            % Energy dissipated by CH nodes after receiving data from all CM nodes
            if(min_dis>0)
               S2(C2(min_dis_cluster).id).E = S2(C2(min_dis_cluster).id).E - ((ERX + EDA)*PS); 
               PACKETS_TO_CH2(r+1)=n-dead2-cluster2+1; 
            end
            S2(i).min_dis=min_dis;
            S2(i).min_dis_cluster=min_dis_cluster;
           
         end
      end
   end
   %hold on;
   
   % For Residual Energy of each nodes per round
       for i=1:n
           RE2(1,r+1)=S2(i).E;
       end
   
   AllRE=0;
   % For Overall Residual Energy of all nodes per round
   for i=1:n
      AllRE=AllRE + S2(i).E;
   end
   OV_RE2(1,r+1)=AllRE;
   
   rcountCHs2=rcountCHs2+countCHs2;
   PACKETS_TO_BS2(r+1)=packets_TO_BS2;
   
end