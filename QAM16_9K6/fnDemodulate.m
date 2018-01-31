function [DiscBit] = fnDemodulate(DiscSymbolIn,hdemod)

%   ���ܣ����������·��Ŀ����һ֡���ݽ����ŵ����ơ�
%   'DiscSymbolIn'  	����һ֡����
%   'M'                 ����M-QAM
%
%   'DiscBit'           ���������bit����

%% Demodulation
zsym = demodulate(hdemod, DiscSymbolIn); % Demodulate signal using 16-QAM.

%% Symbol-to-Bit Mapping
z = de2bi(zsym, 'left-msb');                % Convert integers to bits.% Undo the bit-to-symbol mapping performed earlier.
DiscBit = reshape(z.', numel( z ), 1);	% Convert z from a matrix to a vector.
