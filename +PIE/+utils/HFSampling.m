function [Nx, Ny] = HFSampling(aperture, inCoordsX, inCoordsY, inCoordsZ, outCoordX, outCoordY, outCoordZ, wavelength_nm)


[sr, sc] = size(inCoordsX);

xres = sc;
yres = sr;


xInMin = min(inCoordsX(:));
yInMin = min(inCoordsY(:));
xOutMin = min(outCoordX(:));
yOutMin = min(outCoordY(:));

xInMax = max(inCoordsX(:));
yInMax = max(inCoordsY(:));
xOutMax = max(outCoordX(:));
yOutMax = max(outCoordY(:));

Dix = xInMax - xInMin;
Diy = yInMax - yInMin;

Dox = xOutMax - xOutMin;
Doy = yOutMax - yOutMin;


z_min = min(outCoordZ(:));

lam_um = wavelength_nm/1000;


Nx = Dix*Dox/(z_min*lam_um);
Ny = Diy*Doy/(z_min*lam_um);

if xres < Nx
    fprintf('X Sampling %d is less than the minimum required sampling of %d\n', floor(xres), ceil(Nx));
end

if yres < Ny
     fprintf('Y Sampling %d is less than the minimum required sampling of %d\n', floor(yres), ceil(Ny));
end