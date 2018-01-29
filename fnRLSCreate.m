function [ hRLS ] = fnRLSCreate( TAPs, lamda )
    hRLS.lamda = lamda;
    hRLS.P_last = eye(TAPs)/1;
end
