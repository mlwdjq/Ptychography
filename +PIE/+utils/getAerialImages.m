function [aerialImages,Es] = getAerialImages(E,NAo,Lo_um,NA,lambda_um,z_um,dc_um,df_um,N,offsetAngle)
do_um = lambda_um*z_um/N/dc_um;
k0=2*pi/lambda_um; 
cutoff = NA*k0;
offset = sin(-offsetAngle*NA/NAo/180*pi)*k0;
kmax=pi/dc_um;
[kxm,kym] =meshgrid(linspace(-kmax,kmax,N));
kzm = sqrt(k0^2-(kxm/NA*NAo).^2-(kym/NA*NAo).^2);
CTF = (((kxm-offset).^2+kym.^2)<cutoff^2);
defocus_pha = exp(1i.*df_um.*real(kzm)).*exp(-abs(df_um).*abs(imag(kzm)));
pupil = CTF.*defocus_pha;
spectrum =  PIE.utils.Propagate (E,'fourier',do_um,lambda_um,-1);
spectrum = crop2(spectrum,N,N).*pupil;%imagesc(abs(spectrum));
Es =  PIE.utils.Propagate (spectrum,'fourier',do_um,lambda_um,1);
aerialImages = abs(Es).^2;



%% this approach only work for low-NA
% do_um = lambda_um*z_um/N/dc_um;
% [xp_um,yp_um] = meshgrid(linspace(-do_um*N/2,do_um*N/2,N));
% Rprobe_um = z_um*tan(asin(NA));
% pupil_amp = zeros(N);
% pupil_amp(xp_um.^2+yp_um.^2<= Rprobe_um.^2) = 1;
% kr = sin(atan(sqrt(xp_um.^2+yp_um.^2)/Lo_um));
% phi = atan2(yp_um,xp_um);
% kx = kr.*cos(phi);
% ky = kr.*sin(phi);
% defocus_pha = df_um*2*pi/lambda_um*sqrt(1-kx.^2-ky.^2);
% pupil = pupil_amp.*exp(1i*defocus_pha);
% spectrum =  PIE.utils.Propagate (E,'fourier',do_um,lambda_um,-1);
% spectrum = spectrum.*pupil;
% Es =  PIE.utils.Propagate (spectrum,'fourier',do_um,lambda_um,1);
% aerialImages = abs(Es).^2;