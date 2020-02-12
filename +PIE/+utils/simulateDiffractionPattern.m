function [sqrtInt,Em,measurements] = simulateDiffractionPattern(probe,object,segs,N,propagator,Rpix,H,preShift)
%% simulate diffracted patterns

reconBox = object(Rpix(1)+[1:N],Rpix(2)+[1:N]);
exitWave = reconBox.*probe;
Em = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
sqrtInt = single(abs(Em));
measurements = sqrtInt.^2;


% segment
if ~isempty(segs)
    img = zeros(length(segs{1}));
    out=zeros(1,length(segs));
    for k=1:length(segs)
        seg = logical(segs{k});
        out(k) = sum(sum(measurements(seg)));
        img = img + segs{k}*out(k);
    end
    measurements=img;
    sqrtInt = sqrt(measurements);
end