function [ C ] = Het_Expand_Ring( S, R )
%This function determines the method oh how the ring expands
    C = [];
    New = [];
    for i = 1:size(R,2)
        Tran = Het_Iden(S, R(i));
        for j = 1:size(Tran,2)
            New(size(New,2) + 1) = Tran(j);
        end
        for j = 1:(size(New,2) - 1)
            Tran = path(S, New(j), New(j+1));
            for k = 1:size(Tran,2)
                C(size(C,2) + 1) = Tran(k);
            end
        end
    end
    
    Tran = path(S, C(1), C(size(C,2)));
    for k = 1:size(Tran,2)
        C(size(C,2) + 1) = Tran(k);
    end
    C = unique(C, 'stable');
    for i = 1:size(C,2)-1
        maxd = S(C(i)).ang;
        for j = i+1:size(C,2)
            if(S(C(j)).ang<=maxd)
                maxd = S(C(j)).ang;
                C(j) = C(i) + C(j);
                C(i) = C(j) - C(i);
                C(j) = C(j) - C(i);
            end
        end
    end
    C(size(C,2) + 1) = C(1);   
end