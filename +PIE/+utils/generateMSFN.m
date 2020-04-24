function dMSFN_Phase=generateMSFN(MSFN,nMSFN,minMSFN,maxMSFN)
%  minMSFN minimum frequency
%  maxMSFN maximum frequency
%  nMSFN number of MSFN sampling
if MSFN~=0
    dMSFN_FFT=randn(nMSFN);
    [fx,fy]=meshgrid(1:nMSFN);
    fx=fx-round((nMSFN+1)/2);
    fy=fy-round((nMSFN+1)/2);
    dMSFN_FFT(fx<0|fx.^2+fy.^2<minMSFN^2|fx.^2+fy.^2>maxMSFN^2)=0;
    dMSFN_FFT(fx==0&fy<0)=0;
    dMSFN_FFT=dMSFN_FFT./sqrt(fx.^2+fy.^2);
    dMSFN_FFT(isnan(dMSFN_FFT))=0;
    dMSFN_FFT(isinf(dMSFN_FFT))=0;
    dMSFN_Phase=real(ifft2(ifftshift(dMSFN_FFT)));
    dMSFN_Phase=MSFN*dMSFN_Phase/std(dMSFN_Phase(:));
else
    dMSFN_Phase=0;
end