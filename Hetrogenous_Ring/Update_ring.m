function [ S ] = Update_ring( S, C)
%This function updates the type of the rings.
n = 300;
for i = 1:n
    for j = 1:size(C,2)
        if (i == j)
            S(i).type='ring';
        else
            S(i).type='reg';
        end
    end
end