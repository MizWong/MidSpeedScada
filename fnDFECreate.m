function [ hDFE ] = fnDFECreate( FFTAPs, FBTAPs, RLSLamda, hMOD, hDEM )

hDFE.FFTAPs = FFTAPs;
hDFE.FBTAPs = FBTAPs;
hDFE.RLSLamda = RLSLamda;

hDFE.hRLS = fnRLSCreate( FFTAPs + FBTAPs, RLSLamda );
hDFE.hMOD = hMOD;
hDFE.hDEM = hDEM;

hDFE.wFF = zeros(FFTAPs, 1);
hDFE.wFB = zeros(FBTAPs, 1);

hDFE.e = [];

end
