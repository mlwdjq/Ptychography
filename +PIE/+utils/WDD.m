% reference: Li, P., Edo, T. B., & Rodenburg, J. M. (2014).
% Ptychographic inversion via Wigner distribution deconvolution:
% Noise suppression and probe design. Ultramicroscopy, 147, 106-113.
function [dObjectRecon,dProbeRecon] = WDD(dObjectRecon,dProbeRecon,sqrtInt,N,Rpix,scanSteps,epsilon)
% I_set(u,R), G_set(u,U), H_set(r,U), L_set(r,R)
I_set_u_R = reshape(sqrtInt.^2,N,N,scanSteps,scanSteps);
I_set_R_u= permute(I_set_u_R,[3, 4, 1, 2]);
G_set_u_U = single(zeros(N,N,scanSteps,scanSteps));
G_set_U_u = single(zeros(scanSteps,scanSteps,N,N));
H_set_r_U = G_set_u_U;
H_set_U_r = G_set_U_u;
L_set_r_R = G_set_u_U;
L_set_R_r = G_set_U_u;
O_set_r_R =  G_set_u_U;
P_set_r_R = G_set_u_U;
O_set_u_R = G_set_u_U;
P_set_u_R = G_set_u_U;
O_set_U_u = G_set_U_u;
P_set_U_u = G_set_u_U;
G_set_u_U_s = G_set_u_U;
I_set_R_u_s =I_set_R_u;
D_set_u_U_O =G_set_u_U;
D_set_u_U_P =G_set_u_U;
Pu = G_set_u_U;
Ou = G_set_u_U;
FP = fftshift(fft2(fftshift(dProbeRecon)));
dObjectCrop = single(crop2(dObjectRecon,N,N));
FO = fftshift(fft2(fftshift(dObjectCrop)));
% figure(3),imagesc(abs(dObjectCrop))
[WP,FWP] = PIE.utils.mywigner2(FP,-1);
[WO,FWO] = PIE.utils.mywigner2(FO,1);
%% test
%  H_set_r_U_s = WP.*WO./sqrt(N);
% for m = 1:scanSteps
%     for n = 1:scanSteps
%         G_set_u_U_s(:,:,m,n) = fftshift(fft2(fftshift(H_set_r_U_s(:,:,m,n)))); % from u to r
%     end
% end
% G_set_U_u_s = permute(G_set_u_U_s,[3 4 1 2]);
% for m = 1:N
%     for n = 1:N
%         I_set_R_u_s(:,:,m,n) = ifftshift(ifft2(ifftshift(G_set_U_u_s(:,:,m,n)))); % from R to U
%     end
% end
% I_set_u_R_s= ipermute(I_set_R_u_s,[3, 4, 1, 2]);
% diff = abs(I_set_u_R_s)-I_set_u_R;
% % imagesc(diff(:,:,32,33));colorbar;
% 
% FOs = PIE.utils.deconvolution(FWO,FP,1);
% finish test


% Ps =  pad2(dProbeRecon,length(dObjectRecon),length(dObjectRecon));
% Rpix_flip = max(Rpix(:))-Rpix;
% P_ori_u_U = crop2(fftshift(fft2(fftshift(Ps))),N,N);
% O_ori_u_U = fftshift(fft2(fftshift(dObjectRecon(N/2+1:3*N/2,N/2+1:3*N/2))));
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
%         O_set_r_R(:,:,m,n) = dObjectRecon(Rpix((m-1)*scanSteps+n,1)+[1:N],Rpix((m-1)*scanSteps+n,2)+[1:N]);
%         P_set_r_R(:,:,m,n) = Ps(Rpix_flip((m-1)*scanSteps+n,1)+[1:N],Rpix_flip((m-1)*scanSteps+n,2)+[1:N]);
%         O_set_u_R(:,:,m,n) = fftshift(fft2(fftshift(O_set_r_R(:,:,m,n)))); % from r to u
%         P_set_u_R(:,:,m,n) = fftshift(fft2(fftshift(P_set_r_R(:,:,m,n)))); % from r to u
    end
end
% L_set_R_r = permute(L_set_r_R,[3 4 1 2]);
% H_set_U_r = permute(H_set_r_U,[3 4 1 2]);
% O_set_R_u = permute(O_set_u_R,[3 4 1 2]);
% P_set_R_u = permute(P_set_u_R,[3 4 1 2]);

% for m = 1:N
%     for n = 1:N
%         O_set_U_u(:,:,m,n) = fftshift(fft2(fftshift(O_set_R_u(:,:,m,n)))); % from R to U
%         P_set_U_u(:,:,m,n) = fftshift(fft2(fftshift(P_set_R_u(:,:,m,n)))); % from R to U
%     end
% end
% O_set_u_U = permute(O_set_U_u,[3 4 1 2]);
% P_set_u_U = permute(P_set_U_u,[3 4 1 2]);
% for m = 1:scanSteps
%     for n = 1:scanSteps
%         Xp(:,:,m,n) = ifftshift(ifft2(ifftshift(conj(P_ori_u_U).*P_set_u_U(:,:,m,n))));%% may change the conj
%         Xo(:,:,m,n) = ifftshift(ifft2(ifftshift(conj(O_ori_u_U).*O_set_u_U(:,:,m,n))));%% may change the conj
%     end
% end
Xo = conj(WP).*H_set_r_U./(abs(WP).^2+epsilon)*sqrt(N);
Xp = conj(WO).*H_set_r_U./(abs(WO).^2+epsilon)*sqrt(N);

for m = 1:scanSteps
    for n = 1:scanSteps
        D_set_u_U_O(:,:,m,n) = fftshift(fft2(fftshift(Xo(:,:,m,n))));
        D_set_u_U_P(:,:,m,n) = fftshift(fft2(fftshift(Xp(:,:,m,n))));
    end
end

dObjectRecon = pad2(PIE.utils.deconvolution(D_set_u_U_O,FP,1),length(dObjectRecon),length(dObjectRecon));
dProbeRecon = PIE.utils.deconvolution(D_set_u_U_P,FO,-1);
% figure(3),imagesc(abs(dObjectRecon))

% for m = 1:scanSteps
%     for n = 1:scanSteps
%         Pu(:,:,m,n) = fftshift(fft2(fftshift(Xp(:,:,m,n))));%% may change the conj
%         Ou(:,:,m,n) = fftshift(fft2(fftshift(Xo(:,:,m,n))));%% may change the conj
%     end
% end
% for m = 1:scanSteps
%     for n = 1:scanSteps
%         %         P(m,n) =Pu(round(N/scanSteps*(m-1))+1,round(N/scanSteps*(n-1))+1,m,n);%% may change the conj
%         %         O(m,n) =Ou(round(N/scanSteps*(m-1))+1,round(N/scanSteps*(n-1))+1,m,n);%% may change the conj
%         P(m,n) =Pu(m,n,m,n);%% may change the conj
%         O(m,n) =Ou(m,n,m,n);%% may change the conj
%     end
% end
% dProbeRecon = ifftshift(ifft2(ifftshift(P)));
% dObjectRecon = ifftshift(ifft2(ifftshift(O)));
% q_amp =abs(dProbeRecon);
% q_pha =atan2(imag(dObjectRecon),real(dObjectRecon));
% imagesc(q_amp);colorbar;
% s=1
