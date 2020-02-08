function [q] = WDD(dObjectRecon,dProbeRecon,sqrtInt,N,Rm,Rn,scanSteps)
% I_set(u,R), G_set(u,U), H_set(r,U), L_set(r,R)
I_set_u_R = reshape(sqrtInt.^2,N,N,scanSteps,scanSteps);
I_set_R_u= permute(I_set_u_R,[3, 4, 1, 2]);
G_set_u_U = single(zeros(N,N,scanSteps,scanSteps));
G_set_U_u = single(zeros(scanSteps,scanSteps,N,N));
H_set_r_U = G_set_u_U;
H_set_U_r = G_set_U_u;
L_set_r_R = G_set_u_U;
L_set_R_r = G_set_U_u;
Xa = G_set_u_U;
Q = G_set_u_U;

[Rxpix,Rypix] =meshgrid(Rm,Rn);

for m = 1:N
    for n = 1:N
        G_set_U_u(:,:,m,n) = fftshift(fft2(fftshift(I_set_R_u(:,:,m,n)))); % from R to U
    end
end
G_set_u_U = ipermute(G_set_U_u,[3 4 1 2]);
for m = 1:scanSteps
    for n = 1:scanSteps
        L_set_r_R(:,:,m,n) = ifftshift(ifft2(ifftshift(I_set_u_R(:,:,m,n)))); % from u to r
        H_set_r_U(:,:,m,n) = ifftshift(ifft2(ifftshift(G_set_u_U(:,:,m,n)))); % from u to r
    end
end
L_set_R_r = permute(L_set_r_R,[3 4 1 2]);
H_set_U_r = permute(H_set_r_U,[3 4 1 2]);
for m = 1:scanSteps
    for n = 1:scanSteps
        temp1 = fftshift(fft2(fftshift(dObjectRecon(Rypix(m,n)+[1:N],Rxpix(m,n)+[1:N])))); % from R to U
        temp0 = fftshift(fft2(fftshift(dObjectRecon(Rypix(round((scanSteps+1)/2),...
            round((scanSteps+1)/2))+[-N/2+1:N/2],Rxpix(round((scanSteps+1)/2),round((scanSteps+1)/2))+[-N/2+1:N/2]))));
        Xa(:,:,m,n) = ifftshift(ifft2(ifftshift(conj(temp1).*temp0)));%% may change the conj
    end
end
Xq = conj(Xa).*H_set_r_U./(abs(Xa).^2+eps);
for m = 1:scanSteps
    for n = 1:scanSteps
        Q(:,:,m,n) = fftshift(fft2(fftshift(Xq(:,:,m,n))));%% may change the conj
    end
end
for m = 1:scanSteps
    for n = 1:scanSteps
        Qs(m,n) =Q(round(N/scanSteps*(m-1))+1,round(N/scanSteps*(n-1))+1,m,n);%% may change the conj
    end
end
q = ifftshift(ifft2(ifftshift(Qs)));
q_amp =abs(q);
q_pha =atan2(imag(q),real(q));
imagesc(q_amp);colorbar;
