
function BestSmplPoint = fnFindBestSmplPointTetra(demodSigRRC,N)

rowNum = floor(numel(demodSigRRC)/N);
Rslt = zeros(rowNum,N);
for ii = 1:N
    Rslt(:,ii) =  demodSigRRC(ii : N: ii + N*rowNum - 1);
end

mode = 1;
% varRslt = var(abs(Rslt).^2);
if mode == 1
    sumRslt = sum(abs(Rslt).^2);
    [v,l] = max(sumRslt);
    BestSmplPoint = l;
    
elseif mode == 2
    varRslt = sum(abs(abs(Rslt)-mean(abs(Rslt))));
    [v,l] = min(varRslt);
    BestSmplPoint = l;
end
% plot(sumRslt);grid

end