clear all;
close all;
clc;

%% Parameters
M        =   16;
Rs       =   10e3;
Rb       =   Rs * log2(M);
BB_OSR   =   8;
FsBB     =   Rs * BB_OSR;
IF_OSR   =   36;
FsIF     =   FsBB * IF_OSR;

TFrame   =   0.02;

aRRC     =   0.25;

nTSHead  =   50;
nFTSDvd  =   20;
nFTSLen  =   1;

txcVCO   =   300e6;
txcPPM   =   1;
txcFreq  =   100e3;
% txcFreq  =   txcFreq + (txcVCO / 1e6)*txcPPM;
txcPhi0  =   0;

rxcFreq  =   100e3;
rxcPhi0  =   pi/3;

dSpped   =   350;
fd       =   2;%( txcFreq * dSpped ) / 3e8;
pdb      =   [0 -22.3];
tau      =   [0 5e-6];