%% this script is used for test the ePIE algorithm
% reference: Maiden, A. M., & Rodenburg, J. M. (2009).
% An improved ptychographical phase retrieval algorithm for
% diffractive imaging.?Ultramicroscopy,?109(10), 1256-1262.

addpath('E:\matlab\WZhu_toolbox');
%% parameters setting
lambda_um =13.5e-3; % wavelength
Nobj = 256; % grid number in object plane
object_domian = 200; % object domian, uint: um
Npro = 40; % probe diameter pixel
object_delta =  object_domian/Nobj; % grid size in object plane
dis_obj = 500; % distance from pupil to object plane, unit:um
dis_det = 2000; % distance from object plane to detector plane, unit: um
alpha = 1; % weight factor
beta = 0.1; % weight factor
object_amp = mat2gray(imread('cameraman.tif'))*0.8+0.2; % object amplitude
object_pha=mat2gray(imread('pears.png'));
object_pha= object_pha(1:Nobj,1:Nobj,1);
object_pha = (object_pha-0.5)*2*pi; % object phase
probe_amp = pinhole(Npro,Nobj,Nobj);
probe_pha = zeros(size(object_amp));
Ep = probe_amp.*exp(1i*probe_pha);
Eo = object_amp.*exp(1i*object_pha);
iteration = 1000;
scanSteps = 13; % scanning steps 12*12
Ed = zeros(Nobj,Nobj,scanSteps^2);
scanShift = 8; % scaning shift pixel
[sx,sy] = meshgrid(-floor(scanSteps/2)*scanShift:scanShift:ceil(scanSteps/2)*scanShift-scanShift);
sx(2:2:end,:) = fliplr(sx(2:2:end,:));
sy=sy';
sx=sx';
sxy = [sy(:),sx(:)];
mask = zeros(Nobj);
mask2 = zeros(Nobj);
mask((-floor(scanSteps/2)*scanShift+round(Nobj/2)):(ceil(scanSteps/2)*scanShift-scanShift+round(Nobj/2)),...
    (-floor(scanSteps/2)*scanShift+round(Nobj/2)):(ceil(scanSteps/2)*scanShift-scanShift+round(Nobj/2))) = 1;
mask2((-floor(scanSteps/2)*scanShift+round(Nobj/2)-Npro):(ceil(scanSteps/2)*scanShift-scanShift+round(Nobj/2)+Npro),...
    (-floor(scanSteps/2)*scanShift+round(Nobj/2)-Npro):(ceil(scanSteps/2)*scanShift-scanShift+round(Nobj/2))+Npro) = 1;
Errors = [];
Errors2 = [];
% scanVisualization = object_amp*2; %% scanning visualization
scanVisualization = zeros(Nobj);
overlap = overlapRatio(Npro/2,scanShift); % overlap ratio of two circles
for i = 1:scanSteps^2
    scanVisualization = scanVisualization + circshift(probe_amp,sxy(i,:));
end

simDiff = 1; % simulate diffraction patterns, 1 for yes, 0 for no
%% diffracted pattern simulation using FFT
if simDiff ==1
    Eps=AngSpec(Ep,lambda_um/1000,dis_obj/1000,object_delta/1000); %propagate the probe to the object plane
%          figure(2),imagesc(abs(Eps)),axis equal off tight;
    for i = 1:scanSteps^2
        Ed(:,:,i) = fftshift(fft2(Eo.*circshift(Eps,sxy(i,:))));
        Id = abs(Ed).^2;
        %     figure(2),imagesc(abs(circshift(Eps,[sy(i,j),sx(i,j)]))),axis equal off tight;pause(0.1);
    end
    fprintf('diffraction simulation finished\n');
end
Is = sum(Id,3); % total intensity on detector
%% ePIE
% define initial guesses for probe and object
Eobj = ones(Nobj);
Epro = Ep;
for k = 1:iteration
    % error evaluation
    temp1 = conj(Eobj).*Eo;
    temp2 = abs(Eobj).^2;
    gamma = sum(temp1(mask==1))/sum(temp2(mask==1));
    temp1 = abs(Eo - gamma*Eobj).^2;
    temp2 = abs(Eo).^2;
    Error_sim = sum(temp1(mask==1))/sum(temp2(mask==1));
    Errors(k) = Error_sim;
    temp = 0;
    for j = 1:scanSteps^2
        i=j;
%         if mod(k,2)==0
%             i=scanSteps^2-j+1;
%         end
        phi = Eobj.*circshift(Epro,sxy(i,:)); % exit wave
%                      phi = phi.*mask2;
        Phi = fftshift(fft2(phi));
        Phis = sqrt(Id(:,:,i)).*exp(1i*atan2(imag(Phi),real(Phi))); % updated wave on detector
        temp = temp + abs(sqrt(Id(:,:,i))-abs(Phi)).^2;
        phis = ifft2(ifftshift(Phis)); % updated exit wave
%                      figure(2),imagesc(abs(phis)),axis equal off tight;pause(0.1)
        maxProbes = max(abs(Epro(:)).^2);
        maxObjects = max(abs(Eobj(:)).^2);
        Eobj0=Eobj;
        Eobj = Eobj + alpha*(conj(circshift(Epro,sxy(i,:)))/maxProbes.*(phis-phi));
        
        Epro = Epro + beta*(conj(circshift(Eobj0,-sxy(i,:)*0))/maxObjects.*(phis-phi));
    end
    
    % error evaluation
    Error_exp = sum(temp(mask==1))/sum(Is(mask==1))
    Errors2(k) = Error_exp;
    
    Eobj_amp = abs(Eobj);
    Eobj_pha = atan2(imag(Eobj),real(Eobj));
    Epro_amp = abs(Epro);
    Epro_pha = atan2(imag(Epro),real(Epro));
    Eobj_amp(scanVisualization==0) =NaN;
    Eobj_pha(scanVisualization==0) =NaN;
    Epro_amp(scanVisualization==0) =NaN;
    Epro_pha(scanVisualization==0) =NaN;
    figure(3),subplot(221),imagesc(Eobj_amp),axis equal off tight;title('Sample amplitude');
    figure(3),subplot(222),imagesc(Eobj_pha),axis equal off tight;title('Sample phase');
    figure(3),subplot(223),imagesc(Epro_amp),axis equal off tight;title('Probe amplitude');
    figure(3),subplot(224),imagesc(Epro_pha),axis equal off tight;title('probe phase');drawnow;
    fprintf('%d iterations finished\n',k);
    if k>1&&abs(Errors2(k)-Errors2(k-1))<1e-6
        break;
    end
end


%% plot
% figure(2),imagesc(abs(Eps)),axis equal off tight;hold on;
% figure(2),imagesc(scanVisualization),axis equal off tight;hold on;
% imagesc(probe_amp);
