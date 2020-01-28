function [exitWaveNew,detectorWave] = UpdateExitWave(exitWave,measurement,propagator,dx,wavelength,z)
detectorWave = Propagate(exitWave,propagator,dx,wavelength,z);
correctedWave = sqrt(measurement).*detectorWave./(abs(detectorWave)+eps);
exitWaveNew = Propagate(correctedWave,propagator,dx,wavelength,-z);