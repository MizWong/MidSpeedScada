function [ hDFE ] = fnDFECreate( FFTAPs, FBTAPs, RLSLamda, hMOD, hDEM )

hDFE.FFTAPs = FFTAPs;
hDFE.FBTAPs = FBTAPs;
hDFE.RLSLamda = RLSLamda;

hDFE.hRLS = fnRLSCreate( FFTAPs + FBTAPs, RLSLamda );
hDFE.hMOD = hMOD;
hDFE.hDEM = hDEM;

hDFE.wFF = ones(FFTAPs, 1)/FFTAPs;
hDFE.wFB = ones(FBTAPs, 1)/FBTAPs;

hDFE.e = [];

end
