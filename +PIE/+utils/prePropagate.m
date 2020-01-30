function H = prePropagate (input,propagator,dx,wavelength,z,preShift)
% Pre-propagate a wavefront using a variety of methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: the wavefront to propagate
% propagator: one of 'fourier', 'fresnel' or 'angular spectrum'
% dx: the pixel spacing of the input wavefront
% wavelength: the wavelength of the illumination
% z: the distance to propagate
% H: frequency responds function
% Setup matrices representing reciprocal space coordinates
[ysize,xsize] = size(input);
x = -xsize/2:xsize/2 - 1;
y = -ysize/2:ysize/2 - 1;
fx = x./(dx*xsize);
fy = y./(ysize*dx);
[fx,fy] = meshgrid(fx,fy);
switch propagator
    case 'fourier'
        if z>0
            H = 1;
        else
            H = 0;
        end 
    case 'angular spectrum'
        % Calculate phase distribution for each plane wave component
        w = sqrt(1/wavelength^2 - fx.^2 - fy.^2);
        % exclude evanescent waves
        notEvanescent = imag(w)==0;
        % Compute frequency responds
        H = exp(2i*pi*z*w).*notEvanescent;
    case 'fresnel'
        % Calculate approx phase distribution for each plane wave component
        w = fx.^2 + fy.^2;
        % multiply by phase-shift and inverse transform
        H =exp(-li*pi*z*wavelength*w);
end
if nargin == 6&&preShift==1
    H = ifftshift(H);
end