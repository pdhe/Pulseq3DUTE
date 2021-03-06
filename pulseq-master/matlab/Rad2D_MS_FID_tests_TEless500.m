%% Create a 2D radial FID sequence and export for execution
% 
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
% 
% # Create slice selective RF pulse for imaging.
% # Create readout gradients in both directions 
% # Loop through number of  projections
% # Write the sequence to an open file format suitable for execution on a
% scanner.
% Author: Sairam Geethanath

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assignced default values.
% Pulseq_dir = uigetdir('','Pick the sequences directory');
addpath(genpath('.'));
system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s');

%%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);

%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=256e-3;   %One parameter to test
Nx=256; Ny=256; %resolution along x and y
sliceThickness=5e-3;
TR = 11.5e-3;
TE = 0.5e-3;%minimum 5e-3
dx = fov/Nx;
dy  =dx;

numslices=5; % needed when we do multislice

zstart=sliceThickness*(-numslices/2); %relative to the isocenter 
zend=sliceThickness*(numslices/2);
% radp = get_radkparams(dx,dy,fov);     %?????
%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simultaneously. This
% hierarchical structure means that one event can be used in multiple
% blocks, a common occurrence in MR sequences, particularly in imaging
% sequences. 
%
% First, a slice selective RF pulse (and corresponding slice gradient) can
% be generated using the |makeSincPulse| function.
%
flip=15*pi/180; %flip angle 
[rf, gz] = mr.makeSincPulse(flip,system,'Duration',0.59e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4); %Generate a RF pulse and specify gz

%% Gradients
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{FOV}$$ .  ?????????????????????
% 
% Therefore the area of the readout gradient is $n\Delta k$.
deltak=1/fov;
kWidth = Nx*deltak; % how does kWidth related to gx??

readoutTime = 3.0e-3;
gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime);% Why flatarea is kwidth?

adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);

%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.
preTime=8e-4; %Need to figure this one out later!
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',0.5e-3); % rephase z axis
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*preTime);%???????????????

%% Calculate timing
delayTE=TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf)/2); % this is what we have left for TE and TR after all the operations
% delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
%     - mr.calcDuration(gx) - mr.calcDuration(gzSpoil) - delayTE;
delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
    - mr.calcDuration(gx) - mr.calcDuration(gzSpoil);
delay1 = mr.makeDelay(delayTE);
delay2 = mr.makeDelay(delayTR);

%% Define sequence blocks
% Next, the blocks are put together to form the sequence
Np =32; %Cartesian phase encodes - play around with this number later for evals
Ns = ceil(pi*Np); % Number of spokes required for 360 based on Cartesian Phase encodes
% theta = linspace(0,360,radp.Ns); %This will be replaced with golden angle
dtheta = 360/Ns;
theta = 0:dtheta: (Ns-1)*dtheta; % For debugging
% theta = get_goldenangle(Ns); % For better temporally resolved data acquisition - not meaningful for FID acq

% Np = length(theta); %equivalent to Ns
ktraj = zeros(Ns, adc.numSamples);
scantime=TR*Ns*numslices;
%

%%
for GZstep=zstart:slicethickness:zend
   
    rf.freqOffset=gz.amplitude*GZstep;
    
for np=1:Ns
            seq.addBlock(rf,gz);
            seq.addBlock(gzReph);
            kWidth_projx = kWidth.*cosd(theta(np));
            kWidth_projy = kWidth.*sind(theta(np));

            gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth_projx,'FlatTime',readoutTime);
            gy = mr.makeTrapezoid('y',system,'FlatArea',kWidth_projy,'FlatTime',readoutTime);
            ktraj(np,:) = get_ktraj(gx,gy,adc,1);

            disp(atan2d(gy.flatArea, gx.flatArea));
            seq.addBlock(delay1);
            seq.addBlock(gx,gy,adc);
            seq.addBlock(gzSpoil);
            seq.addBlock(delay2);

end
end

%% Display the sequence 
ktraj = ktraj./max(abs(ktraj(:)))./2;
ktrajs = zeros(size(ktraj,1), size(ktraj,2), 2); 
ktrajs(:,:,1) = real(ktraj);
ktrajs(:,:,2) = imag(ktraj);
figure(1002);
%seq.plot_sg('TimeRange',[0 TR]);
fname = ['Rad2D_MSpei_shortTE_FID', num2str(Ns),'_',num2str(sliceThickness),'ktraj'];
uisave('ktrajs', fname );
%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
fname = [fname, '.seq'];
seq.write(fname)

%%
% % Display the first few lines of the output file
 s=fileread(fname);

 disp(s(1:309))
 scantime=TR*Ns*numslices;
