function output = postPropagate (input,propagator,H,preShift)
% Post-propagate a wavefront using a variety of methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: the wavefront to propagate
% propagator: one of 'fourier', 'fresnel' or 'angular spectrum'
% H: frequency responds function
% output: the propagated wavefront
[ysize,xsize]=size(input);
if nargin==4&&preShift == 1
    switch propagator
        case 'fourier'
            if H>0
                output = fftshift(fft2(fftshift(input)))/sqrt(ysize*xsize);
            else
                output = ifftshift(ifft2(ifftshift(input)))*sqrt(ysize*xsize);
            end
        case 'angular spectrum'
            % Compute FFT of input
            F = fft2(input);
            % multiply FFT by phase-shift and inverse transform
            output = ifft2(F.*H);
        case 'fresnel'
            % Compute FFT
            F = fft2(input);
            % multiply by phase-shift and inverse transform
            output = ifft2(F.*H);
    end
else
    switch propagator
        case 'fourier'
            if H>0
                output = fftshift(fft2(fftshift(input)))/sqrt(ysize*xsize);
            else
                output = ifftshift(ifft2(ifftshift(input)))*sqrt(ysize*xsize);
            end
        case 'angular spectrum'
            % Compute FFT of input
            F = fftshift(fft2(fftshift(input)));
            % multiply FFT by phase-shift and inverse transform
            output = ifftshift(ifft2(ifftshift(F.*H)));
        case 'fresnel'
            % Compute FFT
            F = fftshift(fft2(fftshift(input)));
            % multiply by phase-shift and inverse transform
            output = ifftshift(ifft2(ifftshift(F.*H)));
    end
end