function [sqrtInt,Em,measurements] = simulateDiffractionPattern(probe,object,N,propagator,Rpix,H,preShift)
%% simulate diffracted patterns
reconBox = object(Rpix(1)+[1:N],Rpix(2)+[1:N]);
exitWave = reconBox.*probe;
Em = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
sqrtInt = single(abs(Em));
measurements = sqrtInt.^2;
