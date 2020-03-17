function [sqrtInt,Em,measurements] = simulateDiffractionPattern3D(probe,object,...
    segs,modeNumber,N,propagator,Rpix,H,preShift,do_um,dc_um,lambda_um)
%% simulate diffracted patterns
Em = zeros(N,N,modeNumber);
for m= 1:modeNumber
    reconBox = object(Rpix(1)+[1:N],Rpix(2)+[1:N],m);
    probeTemp = PIE.utils.Propagate(probe(:,:,m),'angular spectrum',do_um,lambda_um,Rpix(3)*1000);
    exitWave = reconBox.*probeTemp;
    Em(:,:,m) = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
%     Em(:,:,m) = PIE.utils.Propagate(EmTemp,'angular spectrum',dc_um,lambda_um,-Rpix(3)*1000);
    sqrtInt = single(abs(Em));
    measurements = sqrtInt.^2;
end
measurements = sum(measurements,3);
sqrtInt = sqrt(measurements);

% segment
if ~isempty(segs)
    img = zeros(length(segs{1}));
    out=zeros(1,length(segs));
    for k=1:length(segs)
        seg = logical(segs{k});
        out(k) = sum(measurements(seg));
        img = img + segs{k}*out(k);
    end
    measurements=img;
    sqrtInt = sqrt(measurements);
end