function [exitWaveNew,detectorWave] = UpdateExitWave(exitWave,sqrtInt,propagator,H,Hm,preShift,dx,wavelength,z)
if nargin <7
    detectorWave = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
    correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps);
    exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm,preShift);
else
    detectorWave = PIE.utils.Propagate(exitWave,propagator,dx,wavelength,z);
    correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps);
    exitWaveNew = PIE.utils.Propagate(correctedWave,propagator,dx,wavelength,-z);
end