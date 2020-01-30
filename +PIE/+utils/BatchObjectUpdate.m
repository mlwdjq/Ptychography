function object = BatchObjectUpdate(exitWaves,probe,object,Rpix,alpha)
[K,L] =size(object);
N = length(probe);
probeSum = zeros(K,L);
waveSum = zeros(K,L);
probeInt = abs(probe).^2;
conjProbe = conj(probe);
J = length(Rpix);
for j = 1:J
    probeSum(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) = probeSum(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) + probeInt;
    waveSum(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) = waveSum(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) + alpha*conjProbe.*exitWaves(:,:,j);
end
object = waveSum./(probeSum+eps);