function [ Anch ] = Het_Mul_Anch( S, no_anch)
%Find one anchor node
    Anch1 = [];
    n = 350;
    for i = 1:n
        d = sqrt((S(i).xd - S(n+2).xd)^2 + (S(i).yd - S(n+2).yd)^2);
        if (d <= (S(i).Power))
            Anch1(size(Anch1,2) + 1) = i;
        end
    end
    for i = 1: size(Anch1,2)
        for j = 1: size(Anch1,2)
            if (S(Anch1(j)).E > S(Anch1(i)).E)
                Anch1(j) = Anch1(i) + Anch1(j);
                Anch1(i) = Anch1(j) - Anch1(i);
                Anch1(j) = Anch1(j) - Anch1(i);
            end
        end
    end
    for i = 1:no_anch
        Anch(i) = Anch1(i);
    end
end