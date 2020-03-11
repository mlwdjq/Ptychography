function [E, amp_near, pha_near,Ex_amp,Ex_phi,Ey_amp,Ey_phi,Ez_amp,Ez_phi,...
    intensity,scatteredOrder_TE_amp,scatteredOrder_TE_pha,scatteredOrder_TM_amp,scatteredOrder_TM_pha]...
    = runSimulator(polarDire,method)
tic,
if nargin<2
    method = 'TEMPESTpr2';
end
switch method
    case 'TEMPESTpr2'
        outh = simulateMaskUsingTEMPESTpr2(1);
    case 'KirchhoffThin'
        outh = simulateMaskUsingStackKirchhoffThick(1);
end
fprintf('near field simulation took %0.1fs\n',toc);

% getDataSeriesText
soh = getDataSeriesHandlesForOutput(outh);
% get data
nP = sqrt(length(getDataSeriesData(soh(1)))); % near field sampling
Ex_amp = reshape(getDataSeriesData(soh(1)),nP,nP);
Ex_phi = reshape(getDataSeriesData(soh(2)),nP,nP);
Ey_amp = reshape(getDataSeriesData(soh(3)),nP,nP);
Ey_phi = reshape(getDataSeriesData(soh(4)),nP,nP);
Ez_amp = reshape(getDataSeriesData(soh(5)),nP,nP);
Ez_phi = reshape(getDataSeriesData(soh(6)),nP,nP);
% E =  Ex_amp.*exp(1i*Ex_phi);
if nargin == 0
    E =  Ex_amp.*exp(1i*Ex_phi)+Ey_amp.*exp(1i*Ey_phi); % average E-field
elseif polarDire == 0
    E =Ex_amp.*exp(1i*Ex_phi);
else
    E =Ey_amp.*exp(1i*Ey_phi);
end
% E =  Ex_amp.*exp(1i*Ex_phi)+Ey_amp.*exp(1i*Ey_phi)+Ez_amp.*exp(1i*Ez_phi); % average E-field
pha_near = atan2(imag(E),real(E)); % near field phase
amp_near = abs(E); % near field amplitude
if nargout>9
    intensity = reshape(getDataSeriesData(soh(7)),nP,nP);
    nS = sqrt(length(getDataSeriesData(soh(8))));
    scatteredOrder_TE_amp = reshape(getDataSeriesData(soh(8)),nS,nS);
    scatteredOrder_TE_pha = reshape(getDataSeriesData(soh(9)),nS,nS);
    scatteredOrder_TM_amp = reshape(getDataSeriesData(soh(10)),nS,nS);
    scatteredOrder_TM_pha = reshape(getDataSeriesData(soh(11)),nS,nS);
end