function [ Tran ] = path( S, first, second)
% This function determines the path to be traversed by the to reach from
% first to second node
n = 350;
Tran = [];
done = 0;
cur = first;
Tran(size(Tran,2) + 1) = first;
while(done~=1)
    indexd = sqrt((S(cur).xd - S(second).xd)^2 + (S(cur).yd - S(second).yd)^2);
    if (indexd <= 80)
        done = 1;
        break
    end
    for j = 1:n
        d1 = sqrt((S(j).xd - S(second).xd)^2 + (S(j).yd - S(second).yd)^2);
        d2 = sqrt((S(j).xd - S(cur).xd)^2 + (S(j).yd - S(cur).yd)^2);
        if (d1<=indexd)&&(d2<=80)
            indexd = d1;
            id = j;
        end
    end
    if (done == 1)
        break
    end
    Tran(size(Tran,2) + 1) = id;
    cur = id;
end
Tran(size(Tran,2) + 1) = second;
end