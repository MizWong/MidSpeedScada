%% 用于scada数传项目
%   要求：
%     1、信道带宽12.5k；
%     2、符号速率9.6k；
%     3、调制方式16qam；
%     4、有效载荷20kbps
%   方案：符号速率9.6k，滚降0.2
%%
clear;
close all;
clc;

%%
SNR_IN = -10 : 2 :10;        %输入信噪比
% SNR_IN = -0;        %输入信噪比
ilen = numel(SNR_IN);
DiffFreq = zeros(1,ilen);
ErrDFreq = zeros(1,ilen);
maxDFreq = zeros(1,ilen);

ErrBit = zeros(1,ilen);
SNR_CHL = zeros(1,ilen);
SNR_Sig = zeros(1,ilen);
DropFrm = zeros(1,ilen);
FrmSyncVar = zeros(1,ilen);
FrmSyncSndVar = zeros(1,ilen);
maxSyncVar =  zeros(1,ilen);
minSyncVar =  zeros(1,ilen)+100;


%% 训练序列
OriTranBit = [1 1 1 1 1 0 0 1 1 0 1 0 1];

%% Setup
M = 16;         % Size of signal constellation
k = log2(M);    % Number of bits per symbol
iStepLen = 12;
n = numel(OriTranBit)*iStepLen*k;        % Number of bits to process
Rs = 9.6e3;
N = 12;
protectLen = 40;
RRC_fs = Rs * N;

ConstPLen = 15;

%%
RRC_115k2Hd = RRC_gen(RRC_fs,N*12,Rs/2,0.2);
RRC_115k2_4k8 = RRC_115k2Hd.numerator;
RRC_115k2_4k8 = RRC_115k2_4k8/sum(RRC_115k2_4k8);
Lp_115k2hd = Lpfilter_window(RRC_fs,5.8e3,6.25e3);
Lp_115k2_5k8_6k25 = Lp_115k2hd.numerator;
Lp_115k2_5k8_6k25 = Lp_115k2_5k8_6k25/sum(Lp_115k2_5k8_6k25);
LP2880hd = LpFlter(2880e3,25e3,55e3);
hd_2880_30_55 = LP2880hd.numerator;
hd_2880_30_55 = hd_2880_30_55/sum(hd_2880_30_55);
Hd = LpFlter(2880e3,10e3,30e3);
hd_fs2880_fpass10_fstop30 = Hd.Numerator;
hd_fs2880_fpass10_fstop30 = hd_fs2880_fpass10_fstop30/sum(hd_fs2880_fpass10_fstop30);
%%
hmod = modem.qammod('M', M, 'SymbolOrder', 'Gray','InputType', 'integer');
hdemod = modem.qamdemod('M', M, 'SymbolOrder', 'Gray');%,'OutputType', 'Bit');

for iTimes = 1:1000
    iTimes
    %% Signal Source
    OriBit = randi([0,1],1,n);     % Random binary data stream
    
    %% Bit-to-Symbol Mapping
    qamMapVal = bi2de(reshape(OriBit, k, length(OriBit) / k).', 'left-msb');% Convert the bits in OriBit into k-bit symbols.
    
    %% Modulation
    qamMapSymbol = modulate(hmod, qamMapVal); % Modulate using 16-QAM.
    
    %% 生成训练序列
    tranSeq = zeros(1,numel(OriTranBit)+1);
    temp = zeros(1,numel(OriTranBit)+1);
    temp(1) = 0;
    for kk = 1:numel(OriTranBit)
        switch(OriTranBit(kk))
            case 0
                temp(kk+1) = mod(temp(kk)-1,4);
            case 1
                temp(kk+1) = mod(temp(kk)+1,4);
        end
    end
    tranSeq(temp==0) = 3+3i;
    tranSeq(temp==1) = -3+3i;
    tranSeq(temp==2) = -3-3i;
    tranSeq(temp==3) = 3-3i;
    
    %% 插入训练序列
    SymbolWithTranSeqLen = numel(qamMapSymbol) / iStepLen * (iStepLen + 1) + 1;    %插入训练序列后总个数
    SymbolWithTranSeq = zeros(1, SymbolWithTranSeqLen);                            %插入训练序列规则：原始数据iStepLen个为一组，将训练序列插入原始序列中，每隔iStepLen个数据插入一个训练序列值
    for kk = 1 : (iStepLen + 1) : (SymbolWithTranSeqLen - 1)
        SymbolWithTranSeq(kk) = tranSeq((kk - 1) / (iStepLen + 1) + 1);
        SymbolWithTranSeq(kk + 1 : kk + iStepLen) = qamMapSymbol((kk - 1) / (iStepLen + 1) * iStepLen + 1 : ((kk - 1) / (iStepLen + 1) + 1) * iStepLen);
    end
    SymbolWithTranSeq(end) = tranSeq(end);
    tranSeqReal = tranSeq;
    
    %% 
    ProtectBit = randi([0,1],protectLen, 1);
    ProtectSymbol = bi2de(reshape(ProtectBit, k, length(ProtectBit) / k).', 'left-msb');% Convert the bits in OriBit into k-bit symbols.
    ProtectMapSymbol = modulate(hmod, ProtectSymbol).'; % Modulate using 16-QAM.
    %     SymbolWithTranSeqPrtct = [ProtectMapSymbol(1:protectLen/4/2).'      SymbolWithTranSeq  ProtectMapSymbol(protectLen/4/2+1:end).'];
    SymbolWithTranSeqPrtct = [ProtectMapSymbol repmat(3+3i,1,ConstPLen) repmat(-3-3i,1,ConstPLen)  repmat(3+3i,1,ConstPLen) SymbolWithTranSeq];%   SymbolWithTranSeq(1:50)];-3-3i 3-3i -3+3i
   
    %% 上采样、成型滤波
    SymbolUp192k = upsample(SymbolWithTranSeqPrtct*N,N);
    SymbolUp192kRRC = conv(SymbolUp192k, RRC_115k2_4k8,'same');
    dataIQ = SymbolUp192kRRC;
    
    %%
    N2 = 2880e3/RRC_fs;
    SymbolUp2880k = upsample(dataIQ*N2,N2);
    SymbolUp2880kLp = conv(SymbolUp2880k,hd_2880_30_55,'same') ;%成型滤波，低通滤波
    
    %%
    fs = 2.88 * 1000000;%采样率2.88M
    fc = 25000; %载波
    t = 0 : 1/fs : (numel(SymbolUp2880kLp) - 1) / fs;
    carryWaveCos = cos(2 * pi * fc * t);
    carryWaveSin = sin(2 * pi * fc * t);
    data_Wave_up = carryWaveCos .* real(SymbolUp2880kLp) - carryWaveSin .* imag(SymbolUp2880kLp);
    
    %% 加入瑞利噪声
    fd = 20;%20、40、114;
    pdb = [0 -22.3];
    tau = [0 5e-6];
    chan = rayleighchan( (1 / fs), fd, tau, pdb);
    data_Ray = filter(chan, data_Wave_up);
    
    seta = (rand(1)*2-1)*2*pi;
    deltF = round((rand(1)*2-1)*400);
    for iSNR = 1:ilen
        
                data_Gauss=data_Ray;
%         data_Gauss = awgn(data_Ray,SNR_IN(iSNR), 'measured'); %加入高斯噪声
        %%
        coefCos = 1;
        coefSin = 1;
        carryWaveCos = coefCos*cos(2 * pi * (fc+deltF) * t + seta);
        carryWaveSin = coefSin*sin(2 * pi * (fc+deltF) * t + seta);
        
        data_disc_cos_1 = data_Gauss .* carryWaveCos;
        data_disc_sin_1 = data_Gauss .* carryWaveSin;
        
        data_disc_IQ = data_disc_cos_1 - i * data_disc_sin_1;
        data_disc = conv(data_disc_IQ,hd_fs2880_fpass10_fstop30,'same') * 2;
        
        %%
        data_Downsingal = conv(data_disc, hd_2880_30_55,'same') ;
        data_Downsample2 = data_Downsingal(1 : N2 : end);%下采样2（15）
        data_DnIQ = conv(data_Downsample2,Lp_115k2_5k8_6k25,'same');
        
        %%
        DemodRRC = conv(data_DnIQ,RRC_115k2_4k8,'same');
        BestSmplPoint = fnFindBestSmplPointTetra(data_DnIQ,N);
        
        %% 位同步
        DemodSymbol = DemodRRC(BestSmplPoint : N : end);
        
%         figure(1);
%         subplot(2,2,1);plot(unwrap(angle(DemodRRC)));title('相位-位同步前')
%         subplot(2,2,3);plot(unwrap(angle(DemodSymbol)));title('相位-位同步后')
%         subplot(2,2,2);plot(abs(diff(unwrap(angle(DemodRRC)))));title('相位差分位同步前')
%         subplot(2,2,4);plot(abs(diff(unwrap(angle(DemodSymbol)))));title('相位差分位同步后')
        
        %% 频差估计 - 位置查找
        PhaseDiff = diff(unwrap(angle(DemodSymbol)));
        PhaseDiffAbs = abs(PhaseDiff);
        ThyPhaseDiffAbs = [zeros(1,ConstPLen-1) pi zeros(1,ConstPLen-1) pi zeros(1,ConstPLen-1)];
        PhaseDiffAbsLen = numel(PhaseDiffAbs);
        ThyPhaseDiffAbsLen = numel(ThyPhaseDiffAbs);
        
        varPhaseDiffAbs = zeros(1,PhaseDiffAbsLen-ThyPhaseDiffAbsLen);
        for kkk = 1:(PhaseDiffAbsLen-ThyPhaseDiffAbsLen)
            varPhaseDiffAbs(kkk) = sum((PhaseDiffAbs(kkk:kkk+ThyPhaseDiffAbsLen-1) - ThyPhaseDiffAbs).^2);
        end
        
        VarLimt = 25;
        constV = pi/2;
        varLoc = [];
        [varVal  varLocAll] = find(varPhaseDiffAbs<VarLimt);
        for jj = 1:numel(varLocAll)
            if varPhaseDiffAbs(varLocAll(jj))<varPhaseDiffAbs(varLocAll(jj)+1) && varPhaseDiffAbs(varLocAll(jj))<varPhaseDiffAbs(varLocAll(jj)-1)
                varLoc = [varLoc varLocAll(jj)];
            end
        end
        if numel(varLoc) == 0
           DropFrm(iSNR)  = DropFrm(iSNR) + 1;
           continue;
        end
        
        varLoc = varLoc(1);
         %% 频差估计 - 相位修正       
        if PhaseDiff(varLoc+ConstPLen-1) > constV
            PhaseDiff(varLoc+ConstPLen-1) = PhaseDiff(varLoc+ConstPLen-1) - pi;
        elseif PhaseDiff(varLoc+ConstPLen-1) < -constV
            PhaseDiff(varLoc+ConstPLen-1) = PhaseDiff(varLoc+ConstPLen-1) + pi;
        end
        if PhaseDiff(varLoc+ConstPLen*2-1) > constV
            PhaseDiff(varLoc+ConstPLen*2-1) = PhaseDiff(varLoc+ConstPLen*2-1) - pi;
        elseif PhaseDiff(varLoc+ConstPLen*2-1) < -constV
            PhaseDiff(varLoc+ConstPLen*2-1) = PhaseDiff(varLoc+ConstPLen*2-1) + pi;
        end
        
        DiffFreq(1,iSNR) = -mean(PhaseDiff(varLoc:varLoc+ConstPLen*3-2))/2/pi*9600;
                
        %% 纠正相位误差
        CrrctSym = DemodSymbol.*exp(j*(DiffFreq(iSNR)/Rs*2*pi*(1:numel(DemodSymbol))));       
%         figure(2);plot(DemodSymbol(10:60),'-*');hold on;plot(CrrctSym(10:60),'-*');hold off;
        figure(2);plot(DemodSymbol,'*');hold on;plot(CrrctSym,'x');hold off;
        %% 帧同步
        [FrmSYNCLoc,val,sndVar] = fnFrmSYNC(CrrctSym,tranSeqReal,iStepLen) ;
        FrmSyncVar(iSNR)  = FrmSyncVar(iSNR) + val;
        FrmSyncSndVar(iSNR)  = FrmSyncSndVar(iSNR) + sndVar;
        
        if val>maxSyncVar(iSNR)
            maxSyncVar(iSNR) = val;
        end
        if sndVar<minSyncVar(iSNR)
            minSyncVar(iSNR) = sndVar;
        end
        
        DemodFrmSym = CrrctSym(FrmSYNCLoc:FrmSYNCLoc+n/4+numel(tranSeqReal)-1);
        
        %% 接收预处理，信道估计
        DemodSymbolCorrected = fnChannelEstimation(DemodFrmSym,tranSeqReal,iStepLen,SymbolWithTranSeqLen);
%         figure(9);plot(DemodFrmSym,'o');hold on;
%         plot(DemodSymbolCorrected,'*');hold off
        
        thyCoef = DemodFrmSym./SymbolWithTranSeq;
        ActCoef = DemodFrmSym./DemodSymbolCorrected;
        figure(10);subplot(2,1,1);plot(abs(thyCoef));hold on;plot(abs(ActCoef));hold off;title('abs')
        subplot(2,1,2);plot(unwrap(angle(thyCoef)));hold on;plot(unwrap(angle(ActCoef)));hold off;title('phase')       
        
        %% Demodulation   Symbol-to-Bit Mapping
        DemodBit = fnDemodulate(DemodSymbolCorrected,hdemod);
        %% 插入训练序列后进行解调，用于误码对比用
        TheorySymbol = demodulate(hdemod, SymbolWithTranSeq);
        
        TheoryBit = de2bi(TheorySymbol, 'left-msb');
        TheoryBit = reshape(TheoryBit.', numel( TheoryBit ), 1);%按先行后列将所有数据变为一列数据，
        
        ErrBit(iSNR) = ErrBit(iSNR) + numel(find(DemodBit-TheoryBit));
    end
   
    ErrDfreqAbs = abs(DiffFreq - deltF);
    ErrDFreq = ErrDFreq + ErrDfreqAbs;
    bigLoc = find((ErrDfreqAbs - maxDFreq)>0);
    maxDFreq(bigLoc) = ErrDfreqAbs(bigLoc);
    avgErrDFreq = ErrDFreq./(iTimes-DropFrm)
    maxDFreq
    Ber = ErrBit/numel(TheoryBit)./(iTimes-DropFrm)
%     [mean(DiffFreq) max(DiffFreq) min(DiffFreq)]
end

ber = ErrBit/numel(TheoryBit)/iTimes;
figure;
semilogy(SNR_IN+20.6,ber)

%% 20180130 16QAM 高斯 SNR_IN = -10:2:10; ConstPLen=15  iStepLen = 12;
avgErrDFreq =   [ 6.0456    4.8814    3.8541    3.0502    2.4216    2.0092    1.6372    1.4016    1.1940    1.0793    0.9291];
maxDFreq =      [29.1002   22.5089   14.2078   12.8556   10.0409   11.5076    7.0236    5.8101    4.6149    4.3356    3.4225];
Ber =           [ 0.0609    0.0326    0.0136    0.0044    0.0010    0.0002    0.0000    0.0000         0         0         0];

%% 20180130 16QAM 瑞利 SNR_IN = -10:2:10; ConstPLen=15  iStepLen = 12;pdb = [0 -22.3];tau = [0 5e-6];
avgErrDFreq_fd2 =    [   6.5664    4.9171    4.2973    3.3499    2.9548    2.4471    2.0911    1.9154    1.8218    1.6960    1.6746];
maxDFreq_fd2 =       [  82.4317   74.0418   53.0288   30.2365   20.9002   29.1789   32.8669   23.4827   28.1130   30.5965   28.8294];
Ber_fd2 =            [   0.0625    0.0338    0.0145    0.0049    0.0013    0.0004    0.0002    0.0001    0.0001    0.0000    0.0000];
avgErrDFreq_fd20 =   [  27.5680   20.8412   18.2779   16.2566   15.6292   14.7696   14.4284   14.3236   13.6310   13.6111   13.5945];
maxDFreq_fd20 =      [ 634.6855  443.8520  554.6643  458.3029  455.7581  450.2868  440.7641  453.8410  227.0293  116.9046  114.6509];
Ber_fd20 =           [   0.0934    0.0632    0.0397    0.0231    0.0154    0.0085    0.0063    0.0038    0.0032    0.0020    0.0019];
DropFrm =            [    58    36    19    10     8     4     3     2     0     0     0];

%% 16QAM  rolloff = 0.2  Lpfilter_window(192e3,5.8e3,6.25e3) 理论同步
SNR_IN = (-10 : 2 :10)+20.6;        %输入信噪比
BER_step11 = [ 0.0475    0.0220    0.0069    0.0014    0.0001    0.0000         0         0         0         0         0];
BER_step7 =  [ 0.0476    0.0216    0.0071    0.0013    0.0001    0.0000         0         0         0         0         0];
BER_step3 =  [ 0.0419    0.0199    0.0066    0.0013    0.0002         0         0         0         0         0         0];
SNR_IN = (-10 : 2 :12)+20.6;        %输入信噪比
BER_stepNone = [ 0.0961    0.0800    0.0654    0.0526    0.0385    0.0255    0.0145    0.0069    0.0025    0.0008    0.0001    6.4e-6];
