function [ delta_w, hRLS ] = fnRLS( hRLS, vu, e, reset )

M = length(vu);

I = eye(M);                 

if reset == 1
    hRLS.P_last = I/1;
end

K = (hRLS.P_last * vu)/(hRLS.lamda + vu'* hRLS.P_last * vu);  
delta_w = conj(K) * e;                 
P = (I - K * vu')* hRLS.P_last/hRLS.lamda;       

hRLS.P_last = P;

end
