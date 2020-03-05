%%
n =33;
Is = permute(I_set_u_R,[3, 2, 1, 4]);
imagesc(Is(:,:,n,n))
Gs = permute(G_set_u_U,[3, 2, 1, 4]);
imagesc(abs(G_set_u_U(:,:,45,33))),axis equal tight
  imagesc(log(abs(Gs(:,:,n,n)))),axis equal tight
%%
Ws = permute(WO,[3, 2, 1, 4]);
imagesc(sqrt(abs(Ws(:,:,31,33))))
imagesc(atan2(imag(Ws(:,:,29,33)),real(Ws(:,:,33,33))))
%%
s =abs(I_set_u_R);
imagesc(s(:,:,21,31));colorbar;
%%
s1 =abs(H_set_r_U_s);
s2 =abs(H_set_r_U);
s =s1-s2;
imagesc(s(:,:,33,33));colorbar;