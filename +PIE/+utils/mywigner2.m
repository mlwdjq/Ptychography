function [w,W] = mywigner2(Exy,dire)
%MYWIGNER: Calculates the Wigner distribution from a column vector
%
%	W  = mywigner(Ex)
%
%	W  = output Wigner distribution
%	Ex = Input electric field (MUST be a column vector
%
%	Notes:
%		W = Int(-inf..inf){E(x+y)E(x-y)exp[2ixy]}
%
%		E(x+y) & E(x-y) are calculated via a FFT (fast Fourier transform) using the
%		shift theorem. The integration is performed via a FFT. Thus it is important
%		for the data to satisfy the sampling theorem:
%		dy = 2*pi/X			X = span of all x-values	dy = y resolution
%		dx = 2*pi/Y			Y = span of all y-values	dx = x resolution
%		The data must be completely contained within the range x(0)..x(N-1) &
%		y(0)..y(N-1) (i.e. the function must fall to zero within this range).
%
%	v1.0
%
%	Currently waiting for update:
%		Remove the fft/ifft by performing this inside the last function calls
%		Allow an arbitrary output resolution
%		Allow an input vector for x (and possibly y).

N = length(Exy);
W = zeros(N,N,N,N);%   Get length of vector
w = zeros(N,N,N,N);%   Get length of vector
[x,y] = meshgrid(ifftshift(((0:N-1)-N/2)*2*pi/(N)));							%   Generate linear vector
[X,Y] = meshgrid((0:N-1)-N/2);
for m =1:N
    for n = 1:N
        EXY1 = fft2( ifft2(Exy));			%   f(u)
        EXY2 = fft2( ifft2(Exy).*exp( dire*1i*x*X(n,m) ).*exp( dire*1i*y*Y(n,m) ));			%   f(u+U)
        W(:,:,m,n) = EXY1.*conj(EXY2);		%   Wigner function
        w(:,:,m,n) = ifftshift(ifft2(ifftshift(EXY1.*conj(EXY2))));		%   Wigner function
%          figure(2),imagesc(abs(EXY1));drawnow;
    end
end
