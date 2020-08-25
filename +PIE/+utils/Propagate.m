function output = Propagate (input,propagator,dx,wavelength,z)
% Propagate a wavefront using a variety of methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: the wavefront to propagate
% propagator: one of 'fourier', 'fresnel' or 'angular spectrum'
% dx: the pixel spacing of the input wavefront
% wavelength: the wavelength of the illumination
% z: the distance to propagate
% output: the propagated wavefront
% Setup matrices representing reciprocal space coordinates
[ysize,xsize] = size(input);
x = -xsize/2:xsize/2 - 1;
y = -ysize/2:ysize/2 - 1;
% x = linspace(-xsize/2,xsize/2,xsize);
% y = linspace(-ysize/2,ysize/2,ysize);
fx = x./(dx*xsize);
fy = y./(ysize*dx);
[fx,fy] = meshgrid(fx,fy);
switch propagator
    case 'fourier'
        if z>0
            output = fftshift(fft2(fftshift(input)))/sqrt(ysize*xsize);
        else
            output = ifftshift(ifft2(ifftshift(input)))*sqrt(ysize*xsize);
        end
    case 'angular spectrum'
        % Calculate phase distribution for each plane wave component
        w = sqrt(1/wavelength^2 - fx.^2 - fy.^2);
        % exclude evanescent waves
        notEvanescent = imag(w)==0;
        % Compute FFT of input
        F = fftshift(fft2(fftshift(input)));
        s = exp(2i*pi*z*w).*notEvanescent;
        s(isnan(s)) = 0;
        % multiply FFT by phase-shift and inverse transform
        output = ifftshift(ifft2(ifftshift(F.*s)));
    case 'fresnel'
        % Calculate approx phase distribution for each plane wave component
        w = fx.^2 + fy.^2;
        % Compute FFT
        F = fftshift(fft2(fftshift(input)));
        % multiply by phase-shift and inverse transform
        output = ifftshift(ifft2(ifftshift(F.*exp(-1i*pi*z*wavelength*w))));
end