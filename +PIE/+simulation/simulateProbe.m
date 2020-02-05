%% this script is used to simulate and save probe 

N= 256;

%% vaccum

probe_amp = pinhole(N/2,N,N);
probe_phase = zeros(N);



%% save probe
probe = probe_amp.*exp(1i*probe_phase);
save('../../data/probe/phaseGrating.mat','probe');