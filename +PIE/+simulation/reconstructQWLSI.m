%% this script is used to reconstruct the phase from QWLSI images

% load aerialImages
Nshift = length(aerialImages);
N = length(aerialImages{1});
Ix = zeros(N,N,Nshift);
Iy =Ix;
Ixs = zeros(Nshift,N^2);
Iys = Ixs;
for i =1:Nshift/2
    Ix(:,:,i) = aerialImages{i};
    Iy(:,:,i) = aerialImages{i+Nshift/2};
    Ixs(i,:) = aerialImages{i}(:);
    Iys(i,:) = aerialImages{i+Nshift/2}(:);
end
%% phase extraction

method = 'random' ;

switch method
    case 'random'
        delta = 0:2*pi/Nshift:2*pi-2*pi/Nshift;
        [dWxs,deltaX] = PIE.utils.RandomShift(Ixs,delta,50);
        [dWys,deltaY] = PIE.utils.RandomShift(Iys,delta,50);
        dWxUnwrapped = reshape(dWxs,N,N);
        dWyUnwrapped = reshape(dWys,N,N);
    case 'PSI'
        dWxUnwrapped = atan2(Ix(:,:,2)-Ix(:,:,4),Ix(:,:,3)-Ix(:,:,1));
        dWyUnwrapped = atan2(Iy(:,:,2)-Iy(:,:,4),Iy(:,:,3)-Iy(:,:,1));
end
% phase unwrap
dWx=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(dWxUnwrapped,255*ones(N));
dWy=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(dWyUnwrapped,255*ones(N));

%% phase reconstruction
shearPercentage  = 0.05*10;
dZ =PIE.utils.Retrieve_LP_iteration(dWx/2/pi, dWy/2/pi, shearPercentage, shearPercentage,ones(255));


%% plot
figure(2),imagesc(dZ);

