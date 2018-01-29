%-------------------------------------------------------------------------%
% Mid Speed SCADA PHY Simulation
% 
%
%
% Miz.Wong ( Miz.Wang@Hytera.com )
% Hytera Co., Ltd.
%
% 2018.3
%-------------------------------------------------------------------------%

clear all;
close all;
figure; fP = 0;

%% Parameters
M        =   16;
Rs       =   10e3;
Rb       =   Rs * log2(M);
BB_OSR   =   8;
FsBB     =   Rs * BB_OSR;
IF_OSR   =   16;
FsIF     =   FsBB * IF_OSR;

coefRRC  =   [-0.000593866582089949,-0.000519178823944003,-0.000307699650483779,2.21345461246340e-06,0.000345155278913264,0.000642233109974616,0.000818884587771089,0.000823014597700411,0.000639151389527122,0.000294953981745590,-0.000142051034444981,-0.000577464052384724,-0.000910387433274922,-0.00105690067175509,-0.000971290586782011,-0.000659806118215797,-0.000183064939475338,0.000354382182865520,0.000825900038269464,0.00111177353217541,0.00112904695336971,0.000855121191294227,0.000338871703412072,-0.000304572965002331,-0.000915753993388131,-0.00132595736521047,-0.00139843265630219,-0.00106657566684062,-0.000359630274176309,0.000591650595676181,0.00157313025628947,0.00232787183764440,0.00261253441190733,0.00225853664726686,0.00122453175737747,-0.000372950590826472,-0.00226059256440084,-0.00404522305227608,-0.00528338406635161,-0.00557283649152548,-0.00464925758407685,-0.00246816369564769,0.000747334927470697,0.00450771296350467,0.00811376852963964,0.0107651950380023,0.0117065485295115,0.0103871815329314,0.00660635327708625,0.000614346285288944,-0.00685448762220582,-0.0146307840696393,-0.0212418585314602,-0.0251118609371785,-0.0248029199554365,-0.0192604225200106,-0.00802204011582344,0.00864660740050646,0.0297168881475067,0.0534748844001800,0.0777062709764141,0.0999650790240016,0.117888685302681,0.129512764698294,0.133538735772974,0.129512764698294,0.117888685302681,0.0999650790240016,0.0777062709764141,0.0534748844001800,0.0297168881475067,0.00864660740050646,-0.00802204011582344,-0.0192604225200106,-0.0248029199554365,-0.0251118609371785,-0.0212418585314602,-0.0146307840696393,-0.00685448762220582,0.000614346285288944,0.00660635327708625,0.0103871815329314,0.0117065485295115,0.0107651950380023,0.00811376852963964,0.00450771296350467,0.000747334927470697,-0.00246816369564769,-0.00464925758407685,-0.00557283649152548,-0.00528338406635161,-0.00404522305227608,-0.00226059256440084,-0.000372950590826472,0.00122453175737747,0.00225853664726686,0.00261253441190733,0.00232787183764440,0.00157313025628947,0.000591650595676181,-0.000359630274176309,-0.00106657566684062,-0.00139843265630219,-0.00132595736521047,-0.000915753993388131,-0.000304572965002331,0.000338871703412072,0.000855121191294227,0.00112904695336971,0.00111177353217541,0.000825900038269464,0.000354382182865520,-0.000183064939475338,-0.000659806118215797,-0.000971290586782011,-0.00105690067175509,-0.000910387433274922,-0.000577464052384724,-0.000142051034444981,0.000294953981745590,0.000639151389527122,0.000823014597700411,0.000818884587771089,0.000642233109974616,0.000345155278913264,2.21345461246340e-06,-0.000307699650483779,-0.000519178823944003,-0.000593866582089949]';
coefIF   =   [-9.71486142686126e-07,-2.52490840134498e-06,-5.15956589130591e-06,-9.11741491425847e-06,-1.45282822708411e-05,-2.13301031515112e-05,-2.91832426566632e-05,-3.73880936708308e-05,-4.48175614068815e-05,-4.98777493783981e-05,-5.05107926080682e-05,-4.42530239107803e-05,-2.83592752072866e-05,6.64332722133572e-20,4.34678941523404e-05,0.000104162861457594,0.000183287816716958,0.000280709858096749,0.000394542341837070,0.000520768254541728,0.000652948990710031,0.000782064905174366,0.000896532654314108,0.000982438859516456,0.00102401980154377,0.00100440281204051,0.000906607277492826,0.000714782608855963,0.000415638414458286,-5.93264313849474e-19,-0.000535598035162358,-0.00118738409584600,-0.00194300316685416,-0.00278047152425816,-0.00366754880661372,-0.00456163124451654,-0.00541024854644544,-0.00615221735218726,-0.00671946756575893,-0.00703951618813806,-0.00703851896063100,-0.00664478608606914,-0.00579260762272626,-0.00442619995375132,-0.00250355990162061,1.74469612064407e-18,0.00308886109919605,0.00674486421327900,0.0109260130742187,0.0155662549159433,0.0205764236385396,0.0258463687940126,0.0312482267484417,0.0366407229122659,0.0418743302294040,0.0467970534213863,0.0512605646449239,0.0551263874131541,0.0582718141062410,0.0605952493112974,0.0620206965373137,0.0625011482793272,0.0620206965373137,0.0605952493112974,0.0582718141062410,0.0551263874131541,0.0512605646449239,0.0467970534213863,0.0418743302294040,0.0366407229122659,0.0312482267484417,0.0258463687940126,0.0205764236385396,0.0155662549159433,0.0109260130742187,0.00674486421327900,0.00308886109919605,1.74469612064407e-18,-0.00250355990162061,-0.00442619995375132,-0.00579260762272626,-0.00664478608606914,-0.00703851896063100,-0.00703951618813806,-0.00671946756575893,-0.00615221735218726,-0.00541024854644544,-0.00456163124451654,-0.00366754880661372,-0.00278047152425816,-0.00194300316685416,-0.00118738409584600,-0.000535598035162358,-5.93264313849474e-19,0.000415638414458286,0.000714782608855963,0.000906607277492826,0.00100440281204051,0.00102401980154377,0.000982438859516456,0.000896532654314108,0.000782064905174366,0.000652948990710031,0.000520768254541728,0.000394542341837070,0.000280709858096749,0.000183287816716958,0.000104162861457594,4.34678941523404e-05,6.64332722133572e-20,-2.83592752072866e-05,-4.42530239107803e-05,-5.05107926080682e-05,-4.98777493783981e-05,-4.48175614068815e-05,-3.73880936708308e-05,-2.91832426566632e-05,-2.13301031515112e-05,-1.45282822708411e-05,-9.11741491425847e-06,-5.15956589130591e-06,-2.52490840134498e-06,-9.71486142686126e-07]';

hMod     =   comm.RectangularQAMModulator( M, 'BitInput', true );
hDem     =   comm.RectangularQAMDemodulator( M, 'BitOutput', true );

T        =   1;

tBB      =   0 : 1/FsBB : T - 1/FsBB;
tIF      =   0 : 1/FsIF : T - 1/FsIF;

txcPPM   =   1;
txcFreq  =   100e3;
txcFreq  =   txcFreq + (txcFreq / 1e6)*txcPPM;
txcPhi0  =   0;

rxcFreq  =   100e3;
rxcPhi0  =   pi/3;

dSpped   =   350;
fd       =   ( txcFreq * dSpped ) / 3e8;
pdb      =   [0 -22.3];
tau      =   [0 5e-6];
RlChan   =   rayleighchan( (1 / FsIF), fd, tau, pdb);

hDFE     =   dfe( 10, 10, rls(0.995), constellation(hMod).' );

hEVM     =   comm.EVM();

%% Transmitter
txBit    =   randi( [0,1], T*Rb, 1 );
txSym    =   step(hMod, txBit );
txSymUp  =   zeros( length(txSym)*BB_OSR, 1 );
txSymUp(1 : BB_OSR : end) = txSym;
txSymUp  =   txSymUp * BB_OSR;
txSymUp  =   FilterConv( txSymUp', coefRRC' )' ;

txSymIF  =   zeros( length(txSymUp)*IF_OSR, 1 );
txSymIF(1 : IF_OSR : end) = txSymUp;
txSymIF  =   txSymIF * IF_OSR;
txSymIF  =   FilterConv( txSymIF', coefIF' )' ;
txSymIF  =   txSymIF .* exp( 1j*( 2*pi*txcFreq*tIF' + txcPhi0 ) );

%% Channel
rxSymIF  =   txSymIF;
rxSymIF  =   filter( RlChan, rxSymIF );
rxSymIF  =   awgn( rxSymIF, 0, 'measured' );

%% Receiver
rxSymIF  =   rxSymIF ./ exp( 1j*( 2*pi*rxcFreq*tIF' + rxcPhi0 ) );
rxSymIF  =   FilterConv( rxSymIF', coefIF' )' ;
rxSymUp  =   zeros( length(rxSymIF) / IF_OSR, 1 );
rxSymUp  =   rxSymIF( 1 : IF_OSR : end );

fP=fP+1;
subplot(3,3,fP);psd(spectrum.periodogram, rxSymUp, 'Fs',FsBB, 'CenterDC',true);title('PSD RX');

fP=fP+1;
subplot(3,3,fP);plot(rxSymUp,'.');title('Pre RRC');

rxSymUp  =   FilterConv( rxSymUp', coefRRC' )';

fP=fP+1;
subplot(3,3,fP);plot(rxSymUp,'.');title('Post RRC');

for ii   =   1:1:BB_OSR
    avePwr(ii) = mean(abs(rxSymUp(ii:BB_OSR:end)));
end
[max T0] =   max(avePwr);
rxSym    =   rxSymUp(T0:BB_OSR:end);

fprintf( 'T0 = %d\r\n', T0 );

fprintf( 'EVM Post-RRC = %f%%\r\n', step(hEVM, txSym, rxSym) );

fP=fP+1;
subplot(3,3,fP);plot(rxSym,'.');title('Post DownSample');

rxSym    =   equalize( hDFE, rxSym, txSym(1:50) );

fP=fP+1;
subplot(3,3,fP);plot(rxSym,'.');title('Post Equalize');

fprintf( 'EVM Post-EQ = %f%%\r\n', step(hEVM, txSym, rxSym) );

(sum( step(hDem, rxSym) ~= txBit ) / length(txBit)) * 100
 
%hs = spectrum.periodogram; figure;psd(hs, rxSym, 'Fs',FSIF, 'CenterDC',true)