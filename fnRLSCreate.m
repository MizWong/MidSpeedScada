function [ hRLS ] = fnRLSCreate( nTAPs, lamda )
% Create RLS system object
% 
% nTaps : Taps number
% lamda : Forgetting factor
%
% Miz.Wong, 2018

    hRLS.lamda = lamda;
    hRLS.P_last = .1*eye(nTAPs);
end
