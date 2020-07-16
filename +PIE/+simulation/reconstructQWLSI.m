%% this script is used to reconstruct the phase from QWLSI images

% load aerialImages
method = 'random' ;
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
%% phase extraction

switch method
    case 'random'
        delta = 0:2*pi/Nshift*2:2*pi-2*pi/Nshift*2;
        [dWxs,deltaX] = PIE.utils.RandomShift(Ixs,delta,50);
        [dWys,deltaY] = PIE.utils.RandomShift(Iys,delta,50);
        deltaX = deltaX/pi*180;
        deltaY = deltaY/pi*180;
        dWxUnwrapped = reshape(dWxs,N,N);
        dWyUnwrapped = reshape(dWys,N,N);
    case 'PSI'
        dWxUnwrapped = atan2(Ix(:,:,2)-Ix(:,:,4),Ix(:,:,3)-Ix(:,:,1));
        dWyUnwrapped = atan2(Iy(:,:,2)-Iy(:,:,4),Iy(:,:,3)-Iy(:,:,1));
        
    case 'fourier'
        xft = (fft(Ix, [], 3));
        yft = (fft(Iy, [], 3));
        dWxUnwrapped =  angle(xft(:,:,2));
        dWyUnwrapped =  angle(yft(:,:,2));
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
scale =4;
sp  = shearPercentage*det0_um/det_um/scale;
dZ =PIE.utils.Retrieve_LP_iteration(dWx/2/pi,dWy/2/pi, sp, sp,ones(N))/scale;
dZs = PIE.utils.DelTilt(dZ);

%% plot
figure(5),imagesc(dWy);colorbar;
figure(2),imagesc(dZs);colorbar;
figure(3),imagesc(Ex_phaN);colorbar;
figure(4),imagesc(dZs-Ex_phaN);colorbar;
%