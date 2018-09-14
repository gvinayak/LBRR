function [ C ] = Het_Iden( S, r )
%This function determines the elements for each node in clockwise and anticlockwise direction
C = [];
Tran1 = [];
Tran2 = [];
n = 350;
k1 = 0;
k2 = 0;
found1 = 0;
found2 = 0;
for i = 1:n
    if(S(i).distance>=S(r).distance)
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
    for i = 1: size(Tran1,2)
        for j = 1: size(Tran1,2)
            if (S(Tran1(j)).E > S(Tran1(i)).E)
                Tran1(j) = Tran1(i) + Tran1(j);
                Tran1(i) = Tran1(j) - Tran1(i);
                Tran1(j) = Tran1(j) - Tran1(i);
            end
        end
    end
    C(1) = Tran1(1);    % Highest Energy element in clockwise direction
end

if(found2==0)
    C(2) = r;
else
    for i = 1: size(Tran2,2)
        for j = 1: size(Tran2,2)
            if (S(Tran2(j)).E > S(Tran2(i)).E)
                Tran2(j) = Tran2(i) + Tran2(j);
                Tran2(i) = Tran2(j) - Tran2(i);
                Tran2(j) = Tran2(j) - Tran2(i);
            end
        end
    end
    C(2) = Tran2(1);    % Highest Energy element in anti-clockwise direction
end
end