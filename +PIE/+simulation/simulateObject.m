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
%% save object
object = object_amp.*exp(1i*object_phase);
save('../../data/object/phaseGrating.mat','object');