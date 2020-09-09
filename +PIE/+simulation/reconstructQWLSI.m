%% this script is used to reconstruct the phase from QWLSI images
 
% load aerialImages
method = 'fourier' ;
if ~strcmp(method,'fourier2D' )
    Nshift = length(aerialImages);
    N = length(aerialImages{1});
    Ix = zeros(N,N,Nshift/2);
    Iy = Ix;
    Ixs = zeros(Nshift/2,N^2);
    Iys = Ixs;
    
    for i =1:Nshift/2
        Ix(:,:,i) = aerialImages{i};
        Iy(:,:,i) = aerialImages{i+Nshift/2};
        Ixs(i,:) = aerialImages{i}(:);
        Iys(i,:) = aerialImages{i+Nshift/2}(:);
    end
end
 
%% try ideal phase shifting
% shifts = linspace(0,2*pi-pi/nInt,nInt/2);
% for i= 1:nInt/2
%     ENx = circshift(EN*exp(1i*shifts(i)),[0,spN])+circshift(EN*exp(-1i*shifts(i)),[0,-spN])+circshift(EN,[spN,0])+circshift(EN,[-spN,0]);
%     Ix(:,:,i) = abs(ENx).^2;
%     Ixs(i,:) =abs(ENx(:)).^2;
%     ENy = circshift(EN,[0,spN])+circshift(EN,[0,-spN])+circshift(EN*exp(1i*shifts(i)),[spN,0])+circshift(EN*exp(-1i*shifts(i)),[-spN,0]);
%     Iy(:,:,i) = abs(ENy).^2;
%     Iys(i,:) =abs(ENx(:)).^2;
% end
%% phase extraction
 
switch method
    case 'random'
        delta = 0:2*pi/Nshift*2:2*pi-2*pi/Nshift*2;
        [dWxs,deltaX] = PIE.utils.RandomShift(Ixs,delta,50);
        [dWys,deltaY] = PIE.utils.RandomShift(Iys,delta,50);
        deltaX = deltaX/pi*180;
        deltaY = deltaY/pi*180;
        unwrapX = unwrap(deltaX/180*pi)/2/pi;
        unwrapX = unwrapX-unwrapX(1)
        dWxUnwrapped = reshape(dWxs,N,N);
        dWyUnwrapped = reshape(dWys,N,N);
    case 'PSI'
        dWxUnwrapped = atan2(Ix(:,:,2)-Ix(:,:,4),Ix(:,:,3)-Ix(:,:,1));
        dWyUnwrapped = atan2(Iy(:,:,2)-Iy(:,:,4),Iy(:,:,3)-Iy(:,:,1));
        
    case 'fourier'
        xft = (fft(Ix, [], 3));
        yft = (fft(Iy, [], 3));
        dWxUnwrapped =  angle(xft(:,:,3));
        dWyUnwrapped =  angle(yft(:,:,3));
    case 'fourier2D'
        Nshift = sqrt(length(aerialImages));
        Ix = zeros(N,N,Nshift);
        Iy = Ix;
        for i =1:Nshift
            for j = 1:Nshift
                Ix(:,:,i) = Ix(:,:,i)+aerialImages{(i-1)*Nshift+j};
                Iy(:,:,j) = Iy(:,:,j)+aerialImages{(i-1)*Nshift+j};
            end
        end
        xft = (fft(Ix, [], 3));
        yft = (fft(Iy, [], 3));
        dWxUnwrapped =  angle(xft(:,:,2));
        dWyUnwrapped =  angle(yft(:,:,2));
end
% phase unwrap
dWx=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(dWxUnwrapped,255*ones(N));
dWy=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(dWyUnwrapped,255*ones(N));
 
%% phase reconstruction
scale = 100;% 2
sp  = shearPercentage*det0_um/det_um/scale;
% spN = sp*N*2;
% dWxs = circshift(Ex_phaN*2*pi,[0,spN])-circshift(Ex_phaN*2*pi,[0,-spN]);
% dWys = circshift(Ex_phaN*2*pi,[spN,0])-circshift(Ex_phaN*2*pi,[-spN,0]);
dZ =PIE.utils.Retrieve_LP_iteration(dWx/2/pi,dWy/2/pi, sp, sp,ones(N))/scale;
% dZ = dZ +aber(:,:,1);
dZs = PIE.utils.DelTilt(dZ);
 
%% plot
xy = linspace(-L_nm/2000,L_nm/2000,N);
figure(5),imagesc(xy,xy,dWy/2/pi);axis tight equal;axis([-inf,inf,-inf,inf,-inf,inf,-0.5,0.5]);
colorbar;xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',16);title('Shearing phase');
figure(2),imagesc(xy,xy,dZs);colorbar;axis tight equal;axis([-inf,inf,-inf,inf,-inf,inf,-0.4,0.1]);
colorbar;xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',16);title('Reconstructed phase');
Ex_phaNs = Ex_phaN;
% Ex_phaNs(Ex_phaNs<-0.15)=NaN;
Ex_phaNs = Ex_phaNs-mean(Ex_phaNs(~isnan(Ex_phaNs)));
figure(3),imagesc(xy,xy,Ex_phaNs);colorbar;axis tight equal;axis([-inf,inf,-inf,inf,-inf,inf,-0.4,0.1]);
colorbar;xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',16);title('Original phase');
residual = dZs-Ex_phaN;
residual =residual -mean(residual(:));
res_crop = residual(80:190,80:190);
 std(res_crop(:))
mask_abs = zeros(N);
mask_sub = zeros(N);
r_abs = 12;
r_abs2 = 25;
r_sub = 60;
shift = 26;
mask_abs(N/2-r_abs+1:N/2+r_abs,N/2-r_abs+1:N/2+r_abs)=1;
mask_sub(N/2-r_sub+1:N/2+r_sub,N/2-r_sub+1:N/2+r_sub)=1;
mask_sub(N/2-r_abs2+1:N/2+r_abs2,N/2-r_abs2+1:N/2+r_abs2)=0;
dW_PS = mean(dZs(mask_sub==1))-mean(dZs(mask_abs==1))
dW0_PS = mean(Ex_phaNs(mask_sub==1))-mean(Ex_phaNs(mask_abs==1))
dWx_PS = (mean(dWx(circshift(mask_abs,[0,-shift])==1))-mean(dWx(circshift(mask_abs,[0,shift])==1)))/2/2/pi;
dWy_PS = (mean(dWy(circshift(mask_abs,[-shift,0])==1))-mean(dWy(circshift(mask_abs,[shift,0])==1)))/2/2/pi;
dWxy_PS = (dWx_PS+dWy_PS)/2;
%    residual(abs(residual)>13*std(residual(:)))=0;
figure(4),imagesc(xy,xy,residual);colorbar;axis tight equal
colorbar;xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',16);title('Residual error');
%
