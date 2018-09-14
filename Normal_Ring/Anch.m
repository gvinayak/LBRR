function [ Tran ] = Anch( S, C )
    %Determines the pattern of transition from anchor node to ring
    A = [];
    n = 300;
    for i = 1:n
        d = sqrt(((S(i).xd - S(n+2).xd)^2) + ((S(i).yd - S(n+2).yd)^2));
        if(d<80)
            A(size(A,2) + 1) = i;
        end
    end
    size(C,2)
    mind = 9999;
    for i = 1:size(A,2)
        for j = 1:size(C,2)
            d = sqrt(((S(A(i)).xd - S(C(j)).xd)^2) + ((S(A(i)).yd - S(C(j)).yd)^2));
            if (d<mind)
                mind = d;
                an.a = i;
                an.r = j;
            end
        end
    end
    for i = 1:size(A,2)
        plot(S(A(i)).xd,S(A(i)).yd,'red .');
    end
    Tran = path(S, A(an.a), an.r);
    for i = 1:size(Tran,2)
        plot(S(Tran(i)).xd,S(Tran(i)).yd,'-');
    end
end