function probe = BatchProbeUpdate(exitWaves,probe,object,Rpix,beta)
N = length(probe);
objectSum = zeros(N);
waveSum = zeros(N);
objectInt = abs(object).^2;
conjObject = conj(object);
J=length(Rpix);
for j=1:J
    objectSum = objectSum + objectInt(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
    waveSum = waveSum + beta*conjObject(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]).*exitWaves(:,:,j);
end
probe = waveSum./(objectSum+eps);