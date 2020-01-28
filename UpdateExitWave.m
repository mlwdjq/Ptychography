function [exitWaveNew,detectorWave] = UpdateExitWave(exitWave,sqrtInt,propagator,H,Hm,preShift,dx,wavelength,z)
if nargin <7
    detectorWave = postPropagate (exitWave,propagator,H,preShift);
    correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps);
    exitWaveNew = postPropagate (correctedWave,propagator,Hm,preShift);
else
    detectorWave = Propagate(exitWave,propagator,dx,wavelength,z);
    correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps);
    exitWaveNew = Propagate(correctedWave,propagator,dx,wavelength,-z);
end