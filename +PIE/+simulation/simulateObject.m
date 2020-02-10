%% this script is used to simulate and save object 

N= 400;

%% vaccum

object_amp = ones(N);
object_phase = zeros(N);

%% amplitude grating
[x,y] = meshgrid(linspace(-1,1,N));

object_amp =4+cos(2*pi*5*x)+cos(2*pi*5*y);
object_amp(object_amp<4) =0;
object_amp(object_amp>=4) =1;

object_phase = zeros(N);
figure,mesh(object_amp)
%% phase grating
[x,y] = meshgrid(linspace(-1,1,N));

object_phase =cos(2*pi*5*x)+cos(2*pi*5*y);

object_amp = ones(N);
figure,mesh(object_phase)

%% Four quadrant phase
object_amp = ones(N);

object_phase = zeros(N);
object_phase(1:N/2,N/2:end)=pi/2;
object_phase(N/2:end,N/2:end)=pi;
object_phase(N/2:end,1:N/2)=pi/2*3;

figure,mesh(object_phase)
%% save object
object = object_amp.*exp(1i*object_phase);
save('../../data/object/fourQuadrantPhase.mat','object');