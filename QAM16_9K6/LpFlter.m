function Hd = LpFlter(Fs,Fpass,Fstop)
%LPFLTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.1 and the DSP System Toolbox 9.3.
% Generated on: 30-Nov-2017 09:53:32

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
% Fs = 48000;  % Sampling Frequency
% 
% Fpass = 2000;            % Passband Frequency
% Fstop = 3000;            % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
