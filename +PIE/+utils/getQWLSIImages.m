function aerialImages = getQWLSIImages(E,NAo,Lo_um,NA,lambda_um,z_um,dc_um,df_um,N,T_um,offsetAngle)
do_um = lambda_um*z_um/N/dc_um;
k0=2*pi/lambda_um; 
cutoff = NA*k0;
offset = sin(-offsetAngle*NA/NAo/180*pi)*k0;
kmax=pi/dc_um;
[kxm,kym] =meshgrid(linspace(-kmax,kmax,N));
kzm = sqrt(k0^2-(kxm/NA*NAo).^2-(kym/NA*NAo).^2);
CTF = (((kxm-offset).^2+kym.^2)<cutoff^2);
defocus_pha = exp(1i.*df_um.*real(kzm)).*exp(-abs(df_um).*abs(imag(kzm)));
s = lambda_um/T_um*z_um; 
QWLSI = exp(1i.*kxm*s)+exp(-1i.*kxm*s)+exp(1i.*kym*s)+exp(-1i.*kym*s);
pupil = CTF.*defocus_pha.*QWLSI;
spectrum =  PIE.utils.Propagate (E,'fourier',do_um,lambda_um,-1);
spectrum = crop2(spectrum,N,N).*pupil;%imagesc(abs(spectrum));
Es =  PIE.utils.Propagate (spectrum,'fourier',do_um,lambda_um,1);
aerialImages = abs(Es).^2;

