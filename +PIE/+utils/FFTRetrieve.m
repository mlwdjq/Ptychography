% Reconstructs a wavefront from its sheared wavefronts using the Fourier
% Transform method. Supports variable shear
% 
% @param {double} dPhiX - Sheared X wavefront
% @param {double} dPhiY - Sheared Y wavefront
% @param {double} dSx   - X shear map
% @param {double} dSy   - Y shear map
% @param {double} dMask - Wavefront mask
% 
% @return {double} dZ   - Reconstructed wavefront

function dZ = FFTRetrieve(dPhiX, dPhiY, dSx, dSy, dMask)
[dSr, dSc] = size(dPhiX);

% Scale sheared wavefronts by shear maps to create proper difference
% functions
dScaledDx = dPhiX / dSc./dSx;
dScaledDy = dPhiY / dSr./dSy;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts

  
dXQ1 = (flipud(dScaledDx));
dXQ2 = (rot90(dScaledDx, 2));
dXQ3 = (fliplr(dScaledDx));
dXQ4 = (dScaledDx);

dYQ1 = (flipud(dScaledDy));
dYQ2 = (rot90(dScaledDy, 2));
dYQ3 = (fliplr(dScaledDy));
dYQ4 = (dScaledDy);

dCompositeDx = [-dXQ2, dXQ1; -dXQ3, dXQ4];
dCompositeDy = [-dYQ2, -dYQ1; dYQ3, dYQ4];

% Transform sheared wavefronts and divide by frequency variable
dM = dSr*2;
dN = dSc*2;

[kr, kc] = meshgrid((0:dM-1)/dM, (0:dN-1)/dN);

% Linear phase
dLP = exp(-2i*pi  * (kr * 0.5 + kc * 0.5) );

dFX = fftshift(fft2(dCompositeDx) .* dLP);
dFY = fftshift(fft2(dCompositeDy) .* dLP);


% dDFx = 1/(dSc * 2);
% dDFy = 1/(dSr * 2);
% 
% 

% delta_f = 1/(N*delta_t);
% fIdxX = (-fix(dM/2):1:fix((dM-1)/2))*dDFx;
% fIdxY = (-fix(dN/2):1:fix((dN-1)/2))*dDFy;
% 
% [dU,dV] = meshgrid(fIdxX,fIdxY);


[h2,w2] = size(dCompositeDx);
[dU,dV] = meshgrid(linspace(-1/2,1/2-1/w2,w2),linspace(-1/2,1/2-1/h2,h2));



dConstructed = (dFX+1i*dFY);

dDivisor = ((dU)+1i*(dV))*2i*pi;
dScaledSpectrum = dConstructed./dDivisor;
% 
% dMultiplier = 1./(dU.^2 + dV.^2) .* (dU - 1i*dV) ;
% dScaledSpectrum = dConstructed.*dMultiplier;



dScaledSpectrum(isinf(dScaledSpectrum) | isnan(dScaledSpectrum)) = 0;


% Inverse transform and mask result
dRec = ifft2(ifftshift(dScaledSpectrum)./dLP ); 

dZ = real(dRec);



dZ = dZ(dSr+1:2*dSr,dSc+1:2*dSc);   
dZ = (dZ).*dMask;




