function [ hDFE, vy, vyDcs ] = fnPLLDFE( hDFE, vu, tsPos, txSym )
% DFE Upgrade
% Miz.Wong

vy = zeros(length(vu), 1);
vyDcs = [];
vyd = zeros(length(vu), 1);
hDFE.e = zeros(length(vu), 1);

% -------------------1/Z-----1/Z-----1/Z ... 1/Z--|
%        |        |       |       |  ...  | 
%        |        |       |       |  ...  |
%        |    /   w1      w2      w3      wN
%        |   /    |       |       |       |
%        |  /     -------------------------
%        | /      |           +           |------------------------------------ y
%        |/       -------------------------          |                       |
%       RLS \     |       |       |       |       Decision                   |
%        |    \   wL     wL-1    wL-2     w1         |                       |
%        |      \ |       |       |       |          |                       |
%        |        |       |       |       |          |                       |
%        |        ---1/Z-----1/Z-----1/Z-----1/Z-----/ -Training Sequence    |
%        |                                      |                            |
%        |                                      |                            |
%        |----------------------------------Calc error-----------------------|


kP = 0.1;
kI = 0.01;
kD = 0;

sIntegrate = 0;
sPhase = 0;
sPhsErrEx = 0;

ePhase = zeros(1, length(vu));

for ii = 1 : 1 : length(vu)
    vu(ii) = vu(ii) * exp(1j*sPhase);
    
    % Step 1 . FFF
    if ii < hDFE.FFTAPs + 1
        yFF(ii) = [zeros(1,hDFE.FFTAPs-ii),vu(1:ii).'] * hDFE.wFF;
    else
        yFF(ii) = vu(ii-hDFE.FFTAPs+1:ii).' * hDFE.wFF;
    end
    
    % Step 2 . FBF
    if ii-1 < hDFE.FBTAPs + 1
        yFB(ii) = [zeros(1,hDFE.FBTAPs-(ii-1)),vyd(1:(ii-1)).'] * hDFE.wFB;
    else
        yFB(ii) = vyd((ii-1)-hDFE.FBTAPs+1:(ii-1)).' * hDFE.wFB;
    end
    
    % Step 3 . Calc yout
    vy(ii) = yFF(ii) + yFB(ii);
    
    % Step 4 . Decision
    tmpDcs = step(hDFE.hDEM, vy(ii));
    vyDcs = [vyDcs, tmpDcs']; 
    
    % Step 5 . Calc error
    if ismember(ii, tsPos)
        vyd(ii) = txSym(ii);
    else
        vyd(ii) = step(hDFE.hMOD, tmpDcs);
    end
    
    hDFE.e(ii) = vyd(ii) - vy(ii);
    
    % Step 6 . Adaptive upgrading
    if ii < hDFE.FFTAPs + 1
        FFInput = [zeros(1,hDFE.FFTAPs-ii),vu(1:ii).'].';
    else
        FFInput = vu(ii-hDFE.FFTAPs+1:ii);
    end
    
    if ii-1 < hDFE.FBTAPs + 1
        FBInput = [zeros(1,hDFE.FBTAPs-(ii-1)),vyd(1:(ii-1)).'].';
    else
        FBInput = vyd((ii-1)-hDFE.FBTAPs+1:(ii-1));
    end
    
    [deltaW, hDFE.hRLS] = fnRLS( hDFE.hRLS, [FFInput;FBInput], hDFE.e(ii), 0 );
    
    deltaW_FF = deltaW(1:hDFE.FFTAPs);
    deltaW_FB = deltaW(end-hDFE.FBTAPs+1:end);
    
    % Step 7 . Upgrade weight
    hDFE.wFF = hDFE.wFF + deltaW_FF;
    hDFE.wFB = hDFE.wFB + deltaW_FB;
    
    ePhase(ii) = phase( vyd(ii) / vy(ii) );
    sIntegrate = sIntegrate + ePhase(ii)*kI;
    sPhase = sPhase + kP * ePhase(ii) + sIntegrate + kD * (ePhase(ii) - sPhsErrEx);
end

hDFE.ePhase = ePhase;

end

