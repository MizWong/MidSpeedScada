function [ hDFE ] = fnDFECreate( FFTAPs, FBTAPs, RLSLamda, hDEM )

hDFE.FFTAPs = FFTAPs;
hDFE.FBTAPs = FBTAPs;
hDFE.RLSLamda = RLSLamda;

hDFE.hFFRLS_I = fnRLSCreate( FFTAPs, RLSLamda );
hDFE.hFFRLS_Q = fnRLSCreate( FFTAPs, RLSLamda );
hDFE.hFBRLS_I = fnRLSCreate( FBTAPs, RLSLamda );
hDFE.hFBRLS_Q = fnRLSCreate( FBTAPs, RLSLamda );
hDFE.hDEM = hDEM;

end
