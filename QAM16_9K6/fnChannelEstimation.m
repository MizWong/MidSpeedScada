function [DiscSymbolCorrected] = fnChannelEstimation(DiscSymbol,tranSeqReal,iStepLen,SymbolWithTranSeqLen)

%   ���ܣ����������·��Ŀ����һ֡���ݽ����ŵ����ơ�
%   'DiscSymbol'            ����һ֡����
%   'tranSeqReal'           ����ѵ������
%   'iStepLen'            	����ѵ�����в��벽������������������һ��ѵ������ֵ���򲽳�Ϊ5
%   'SymbolWithTranSeqLen' 	�������ѵ�����к󳤶�
%
%   'DiscSymbolCorrected'   ��������ŵ����ƺ������

RealTranSeqVal = DiscSymbol(1 : (iStepLen + 1) : SymbolWithTranSeqLen);     %ȡ�����յ��������ж�Ӧѵ�����е�ֵ
CorrReal = RealTranSeqVal ./ tranSeqReal;                                   %����ԭֵ���õ�������仯ϵ��

TranSeqLoc = 1 : (iStepLen + 1) : SymbolWithTranSeqLen;
% CorrCorrected = interp1(TranSeqLoc, CorrReal, 1 : SymbolWithTranSeqLen, 'spline');	%��ϵ�������ڲ崦���õ��������ݵ�ϵ���仯ֵ
CorrCorrected = interp1(TranSeqLoc, CorrReal, 1 : SymbolWithTranSeqLen, 'linear');	%��ϵ�������ڲ崦���õ��������ݵ�ϵ���仯ֵ
DiscSymbolCorrected = DiscSymbol ./ CorrCorrected;                                  % ʹ���ڲ�õ���ϵ���仯ֵ�Խ��յ������ݽ�������

% % figure;
% subplot(1,2,1)
% plot(TranSeqLoc,real(CorrReal),'--',TranSeqLoc,real(CorrReal),'+',1 : SymbolWithTranSeqLen,real(CorrCorrected),'x') % �����Լ����η���ʽ�ڲ��ͼ
% legend('ѵ�����е����ϵ��','���Բ�ֵ��õ��������л���ϵ��')
% subplot(1,2,2)
% plot(TranSeqLoc,imag(CorrReal),'--',TranSeqLoc,imag(CorrReal),'+',1 : SymbolWithTranSeqLen,imag(CorrCorrected),'x') % �����Է���ʽ��spline�ڲ��ͼ
% legend('ѵ�����е����ϵ��','���Բ�ֵ��õ��������л���ϵ��')
% hold off
