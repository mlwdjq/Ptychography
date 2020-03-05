% reference: Li, P., Edo, T. B., & Rodenburg, J. M. (2014).
% Ptychographic inversion via Wigner distribution deconvolution:
% Noise suppression and probe design. Ultramicroscopy, 147, 106-113.
function O = deconvolution(D,FP,dire)
N = length(D);
x = conj(D(round((N+1)/2),round((N+1)/2),:,:))./...
sqrt(D(round((N+1)/2),round((N+1)/2),round((N+1)/2),round((N+1)/2)));
xs=permute(x,[3, 4, 1, 2])'; % get FO
xsShift = zeros(N,N,N,N);%   Get length of vector
xsShift2 = zeros(N,N,N,N);%   Get length of vector

[x,y] = meshgrid(ifftshift(((0:N-1)-N/2)*2*pi/(N)));		%   Generate linear vector
[X,Y] = meshgrid((0:N-1)-N/2);
s1 = 0;
s2 = 0;
for m =1:N
    for n = 1:N
        xsShift(:,:,m,n) = fft2( ifft2(xs).*exp( dire*1i*x*X(n,m)/2*2 ).*exp( dire*1i*y*Y(n,m)/2*2 ));			%   f(u+U)
        xsShift2(:,:,m,n) = abs(xsShift(:,:,m,n)).^2;
        if FP(m,n)<0.001
            s1 = s1 + (xsShift(:,:,m,n)).*D(:,:,m,n);
            s2 = s2 + xsShift2(:,:,m,n);
        end
    end
end
FO = s1./s2;
O=(ifftshift(ifft2(ifftshift(FO))));

