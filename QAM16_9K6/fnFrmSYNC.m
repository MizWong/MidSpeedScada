function [FrmBegLoc,val,sndVar] = fnFrmSYNC(DemodSymbol,tranSeqReal,iStepLen)

dataInLen = numel(DemodSymbol);
tranSeqRealLen = numel(tranSeqReal);

RsltLen = dataInLen -  (iStepLen+1) * (tranSeqRealLen-1) ;
Rslt = zeros(RsltLen,tranSeqRealLen);
for ii = 1:RsltLen
    Rslt(ii,:) =  DemodSymbol(ii : iStepLen+1 : ii + (iStepLen+1)*tranSeqRealLen - 1);
end

DiffAngleRslt =  diff(angle(Rslt(:,:)).').';
SmallThanPiRslt = find(DiffAngleRslt < -pi);
DiffAngleRslt(SmallThanPiRslt) = DiffAngleRslt(SmallThanPiRslt) + 2 * pi;
BigThanPiRslt = find(DiffAngleRslt > pi);
DiffAngleRslt(BigThanPiRslt) = DiffAngleRslt(BigThanPiRslt) - 2 * pi;

DiffAngleThy =  diff(angle(tranSeqReal(:,:)).').';
SmallThanPiRslt = find(DiffAngleThy < -pi);
DiffAngleThy(SmallThanPiRslt) = DiffAngleThy(SmallThanPiRslt) + 2 * pi;
BigThanPiRslt = find(DiffAngleThy > pi);
DiffAngleThy(BigThanPiRslt) = DiffAngleThy(BigThanPiRslt) - 2 * pi;


VarR = sum((DiffAngleRslt-DiffAngleThy).^2,2);
% CorrR = sum(DiffAngleRslt.*DiffAngleThy,2);
% plot(VarR);grid

[val,Loc] = min(VarR);
sortVar = sort(VarR);
FrmBegLoc = Loc;
sndVar = sortVar(2);
end