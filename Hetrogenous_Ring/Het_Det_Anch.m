function [ Anch ] = Het_Det_Anch( S, C )
%Detects the anchor nodes for heterogenous network
    fix_rate = 200;
    buffer = 50;
    rec_rate = size(S,2) - 2 - size(C,2);
    
    if(fix_rate >= rec_rate)
        f1 = @Het_One_Anch;
        Anch = f1(S);
        
    else
        rate = rec_rate - fix_rate;
        no_anch = ceil(rate/buffer); 
        f1 = @Het_Mul_Anch;
        Anch = f1(S, no_anch);
    end
end