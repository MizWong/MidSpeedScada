function [DiscSymbolCorrected] = fnChannelEstimation(DiscSymbol,tranSeqReal,iStepLen,SymbolWithTranSeqLen)

%   功能：针对无线链路项目，对一帧数据进行信道估计。
%   'DiscSymbol'            输入一帧数据
%   'tranSeqReal'           输入训练序列
%   'iStepLen'            	输入训练序列插入步长，若六个符号中有一个训练序列值，则步长为5
%   'SymbolWithTranSeqLen' 	输入插入训练序列后长度
%
%   'DiscSymbolCorrected'   输出进行信道估计后的数据

RealTranSeqVal = DiscSymbol(1 : (iStepLen + 1) : SymbolWithTranSeqLen);     %取出接收到的数据中对应训练序列的值
CorrReal = RealTranSeqVal ./ tranSeqReal;                                   %除以原值，得到各个点变化系数

TranSeqLoc = 1 : (iStepLen + 1) : SymbolWithTranSeqLen;
% CorrCorrected = interp1(TranSeqLoc, CorrReal, 1 : SymbolWithTranSeqLen, 'spline');	%对系数进行内插处理，得到所有数据的系数变化值
CorrCorrected = interp1(TranSeqLoc, CorrReal, 1 : SymbolWithTranSeqLen, 'linear');	%对系数进行内插处理，得到所有数据的系数变化值
DiscSymbolCorrected = DiscSymbol ./ CorrCorrected;                                  % 使用内插得到的系数变化值对接收到的数据进行修正

% % figure;
% subplot(1,2,1)
% plot(TranSeqLoc,real(CorrReal),'--',TranSeqLoc,real(CorrReal),'+',1 : SymbolWithTranSeqLen,real(CorrCorrected),'x') % 将线性及三次方程式内插绘图
% legend('训练序列点畸变系数','线性差值后得到整个序列畸变系数')
% subplot(1,2,2)
% plot(TranSeqLoc,imag(CorrReal),'--',TranSeqLoc,imag(CorrReal),'+',1 : SymbolWithTranSeqLen,imag(CorrCorrected),'x') % 将线性方程式及spline内插绘图
% legend('训练序列点畸变系数','线性差值后得到整个序列畸变系数')
% hold off
