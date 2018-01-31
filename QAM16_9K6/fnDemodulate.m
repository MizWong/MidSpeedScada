function [DiscBit] = fnDemodulate(DiscSymbolIn,hdemod)

%   功能：针对无线链路项目，对一帧数据进行信道估计。
%   'DiscSymbolIn'  	输入一帧数据
%   'M'                 输入M-QAM
%
%   'DiscBit'           输出解调后的bit数据

%% Demodulation
zsym = demodulate(hdemod, DiscSymbolIn); % Demodulate signal using 16-QAM.

%% Symbol-to-Bit Mapping
z = de2bi(zsym, 'left-msb');                % Convert integers to bits.% Undo the bit-to-symbol mapping performed earlier.
DiscBit = reshape(z.', numel( z ), 1);	% Convert z from a matrix to a vector.
