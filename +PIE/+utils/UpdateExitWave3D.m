function [exitWaveNew,detectorWave] = UpdateExitWave3D(exitWave,sqrtInt,propagator,H,Hm,preShift,dc_um,lambda_um,Rpix,dx,z)
if nargin <10
    detectorWave = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
    detectorWave = PIE.utils.Propagate(detectorWave,'angular spectrum',dc_um,lambda_um,-Rpix(3)*1000);
    correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps);
    correctedWave = PIE.utils.Propagate(correctedWave,'angular spectrum',dc_um,lambda_um,Rpix(3)*1000);
    exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm,preShift);
else
    detectorWave = PIE.utils.Propagate(exitWave,propagator,dx,lambda_um,z);
    detectorWave = PIE.utils.Propagate(detectorWave,'angular spectrum',dc_um,lambda_um,-Rpix(3)*1000);
    correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps);
    correctedWave = PIE.utils.Propagate(correctedWave,'angular spectrum',dc_um,lambda_um,Rpix(3)*1000);
    exitWaveNew = PIE.utils.Propagate(correctedWave,propagator,dx,lambda_um,-z);
end