function [ hRLS ] = fnRLSCreate( TAPs, lamda )
    hRLS.lamda = lamda;
    hRLS.P_last = .1*eye(TAPs);
end
