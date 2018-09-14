function [ C ] = Iden2( S, r )
%This function determines the elements for each node in clockwise and anticlockwise direction
C = [];
Tran1 = [];
Tran2 = [];
n = 300;
k1 = 0;
k2 = 0;
found1 = 0;
found2 = 0;
for i = 1:n
    if(S(i).distance<=S(r).distance)
        d = sqrt((S(i).xd - S(r).xd)^2 + (S(i).yd - S(r).yd)^2);
        if(d<80)
            %Clockwise Direction
            if(S(i).ang <= S(r).ang)
                found1 = 1;
                k1 = k1 + 1;
                Tran1(k1) = i;
            %Anti clockwise Direction
            else
                found2 = 1;
                k2 = k2 + 1;
                Tran2(k2) = i;
            end
        end
    end
end

if(found1 == 0)
    C(1) = r;
else
    C(1) = Tran1(randi(k1));
end

if(found2==0)
    C(2) = r;
else
    C(2) = Tran2(randi(k2));
end
end