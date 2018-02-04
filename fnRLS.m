function [ delta_w, hRLS ] = fnRLS( hRLS, vu, e, reset )
% Upgrade RLS system object and calculate delta weight
% 
% vu : u
% e  : error
%
% Miz.Wong, 2018

M = length(vu);

I = eye(M);                 

if reset == 1
    hRLS.P_last = .1*I;
end

% Step 1 . Kalman gain
K = (hRLS.P_last * vu)/(hRLS.lamda + vu'* hRLS.P_last * vu);  
% Step 2 . Calc delta weight
delta_w = conj(K) * e;      
% Step 3 . Upgrade P
P = (hRLS.P_last-K*vu'*hRLS.P_last)/hRLS.lamda; 

hRLS.P_last = P;

end
