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

function dZ = CDFTRetrieve(dPhiX, dPhiY, dSx, dSy, dMask)
[dSr, dSc] = size(dPhiX);

% Scale sheared wavefronts by shear maps to create proper difference
% functions
dScaledDx = dPhiX / dSc./dSx;
dScaledDy = dPhiY / dSr./dSy;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts

  
dXQ1 = padTL(flipud(dScaledDx));
dXQ2 = padTL(rot90(dScaledDx, 2));
dXQ3 = padTL(fliplr(dScaledDx));
dXQ4 = padTL(dScaledDx);

dYQ1 = padTL(flipud(dScaledDy));
dYQ2 = padTL(rot90(dScaledDy, 2));
dYQ3 = padTL(fliplr(dScaledDy));
dYQ4 = padTL(dScaledDy);

dCompositeDx = [-dXQ2, dXQ1; -dXQ3, dXQ4];
dCompositeDy = [-dYQ2, -dYQ1; dYQ3, dYQ4];

% dCompositeDx = fourierShift2(dCompositeDx, [.5, .5]);
% dCompositeDy = fourierShift2(dCompositeDy, [.5, .5]);

% Transform sheared wavefronts and divide by frequency variable
dFX = (fft2(dCompositeDx));
dFY = (fft2(dCompositeDy));

dSc = dSc + 1;
dSr = dSr + 1;

dDFx = 1/(dSc * 2);
dDFy = 1/(dSr * 2);

dM = dSr*2;
dN = dSc*2;

% delta_f = 1/(N*delta_t);
fIdxX = (-fix(dM/2):1:fix((dM-1)/2))*dDFx;
fIdxY = (-fix(dN/2):1:fix((dN-1)/2))*dDFy;

[dU,dV] = meshgrid(fIdxX,fIdxY);

dScaledSpectrum = (dFX+1i*dFY)./(fftshift(dU)+1i*fftshift(dV))./2/pi/1i;



dScaledSpectrum(isinf(dScaledSpectrum) | isnan(dScaledSpectrum)) = 0;

% Inverse transform and mask result
dRec = ifft2((dScaledSpectrum)); 
dZ = real(dRec);



dZ = dZ(dSr+1:2*dSr - 1,dSc+1:2*dSc - 1);   
dZ = (dZ).*dMask;

function out = padTL(in)

[sr, sc] = size(in);
out = zeros(sr + 1, sc + 1);
out(2:end, 2:end) = in;


