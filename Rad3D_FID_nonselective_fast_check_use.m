%% Create a 3D radial FID sequence and export for execution
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
% Author: Sairam Geethanath - tricky with the slice sel gradient strength
%         Peidong He
% :)

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assignced default values.
% Pulseq_dir = uigetdir('','Pick the sequences directory');
addpath(genpath('.'));
dw = 10e-6; %dwell time 
rfdt = 10e-6; %RF dead time
rfrt = 10e-6; %RF ringdown time
system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s', 'rfRingdownTime', rfrt, 'rfDeadtime',rfdt);

%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);

%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
% get_radkparams returns the total spokes, including polar and azimuthal, 
% needed for defined unaliased FOV. It follows the projection acquisition 
% theory to calculate the required spokes.

fov=256e-3;   %One parameter to test
Nx=128; Nz = 128;% Will work with lesser for now - change it back to 128
sliceThickness=30e-3;%was 10e-3 before
TR = 20e-3; %repetition time

dx = fov/Nx; %spatial resolution
dz  = sliceThickness/Nz; %spatial resolution at z direction
radp = get_radkparams(dz,dx,fov,'3D'); %Contains number of spokes, Ntheta an Nphi in the polar coordinates
wr_traj =1;


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
flip=15*pi/180; % flip angle= 15degree
st = 30e-3; %slice thickness
% [rf, gz] = mr.makeSincPulse(flip,system,'Duration',1e-3,...
%     'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
% [rf,gz]=mr.makeBlockPulse(flip, 'Duration',20e-6,'SliceThickness',st, 'system', system);
[~,gz] = mr.makeBlockPulse(pi/30,'Duration',dw, 'SliceThickness',st ,'system', system);

% Rectangle RF does not need Gz for slice selection
% Gz  here is for spoiling gradient later


rf_dur = 1*dw; %RF duration
%[rf] = mr.makeBlockPulse(pi/30,'Duration',rf_dur, 'system', system);
[rf] = mr.makeBlockPulse(flip,'Duration',rf_dur, 'system', system); %250khz bandwidth
%% Gradients
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{FOV}$$
% 
% Therefore the area of the readout gradient is $n\Delta k$.

deltak=1/fov; 
kWidth = Nx*deltak;% in radial case, it is kmax

readoutTime = 6.4e-3/5;% sw for this is 128/6.4=20khz which is not good, should be 100khz.
gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime, 'system', system);% sample during flattop

%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.
preTime=8e-4; %Need to figure this one out later!
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',0.5e-3); %rephase
gxSpoil = mr.makeTrapezoid('x',system,'Area',gz.area*2,'Duration',3*preTime);
gySpoil = mr.makeTrapezoid('y',system,'Area',gz.area*2,'Duration',3*preTime);
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',7*preTime); %changed this from 7*preTime
% should duration be the same? not necessarily 
% The code for gradient spoiling is being demonstrated here.
% Gx, Gy and Gz can be applied together 


%% Calculate timing
% delayTE=TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf)/2);
%delayTE=TE - (mr.calcDuration(rf)/2);
delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ... %Tricky here, not consistent
    - mr.calcDuration(gx) - mr.calcDuration(gzSpoil);

delayTRps = 0.8.*delayTR;%to avoid gradient heating
delayTRps = delayTRps - mod(delayTRps, 1e-5);% get rid of the remainder 
delayTRs = 0.2.*delayTR;
delayTRs = delayTRs - mod(delayTRs, 1e-5);


%delay1 = mr.makeDelay(delayTE);
delay2 = mr.makeDelay(delayTRps);
delay3 = mr.makeDelay(delayTRs);

%Delay TR is basically TR subtract by the time needed for RF, readout and spoiling. 
%The delay TR is being splitted to 80% and 20%. 80% of time is between 
%read out and spoiling, 20% of time is placed between spoiling and next RF


%% Changing values of spokes to save time

radp.Ns=16384;
radp.Ntheta=radp.Ns^0.5;
radp.Nphi=radp.Ns^0.5;

% The ideal case required 3748066() spokes which is not clinical applicable. 
% Therefore we decrease the number of spokes here and use under-sample 
% reconstruction later  
 


%% Define spoke angles - polar (theta); azimuthal (phi)

theta = linspace(0,179,radp.Ntheta); %Polar
phi = linspace(0,359,radp.Nphi); %azimuthal

% square root theta and phi
% When get the numbers of spokes, start to distribute them 
% in polar angle and azimuthal angle.

% Pre-registered the memory space for k-trajectory here.

%% Define sequence blocks
ktraj = zeros(radp.Ns, adc.numSamples,3);
np = 0;

%% Create AFP pulse Adiabatic pulses
% beta = 800; %rad/s
% myu = 4.9;
% A0 = 2;
% Tp=gz.flatTime;
% N=Nx;
% t = linspace(-Tp/2, Tp/2,N); %AFP
% t = linspace(-Tp, Tp,N); %AFP
% 
% [B1,A,w1] = make_HSI(A0,beta,myu,t,1); %AFP
% [rf_full] = mr.makeArbitraryRf(B1,flip, 'system', system);%%



%% Actual loop that writes the sequence
for nt=1:radp.Ntheta % 0 - pi
     for ph = 1: radp.Nphi % 0-2pi
            % Excitation

            seq.addBlock(rf);
            % Estimate required kspace extents
            kWidth_projx = kWidth.*sind(theta(nt)).*cosd(phi(ph));
            kWidth_projy = kWidth.*sind(theta(nt)).*sind(phi(ph));
            kWidth_projz = kWidth.*cosd(theta(nt));

            % Determine gradient waveforms for each direction
            gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth_projx,'FlatTime',readoutTime);
            gy = mr.makeTrapezoid('y',system,'FlatArea',kWidth_projy,'FlatTime',readoutTime);
            gz = mr.makeTrapezoid('z',system,'FlatArea',kWidth_projz,'FlatTime',readoutTime);
            
            % Determine the kspace trajectory and store it for
            % reconstruction
            G.gx = gx;G.gy = gy;G.gz = gz;
            np = np+1;
            ktraj(np,:,:) = get_ktraj_v1(G,adc,0);
           
%             seq.addBlock(delay1);
            seq.addBlock(gx,gy,gz,adc);
            seq.addBlock(delay2);
%             seq.addBlock(gzSpoil);%
            seq.addBlock(gxSpoil, gySpoil, gzSpoil);
            delay3 = mr.makeDelay(delayTRs);
            
     end
      disp(nt/radp.Ntheta);
end
%% Normalize the trajectory  and display the sequence
ktrajs = ktraj./max(abs(ktraj(:)))./2; 


%??????????????????????????????????????
% figure();
% seq.plot_sg('TimeRange',[0 TR]);
%________________________________________________Why the plot_sg does not
%work? GPI doesnt take stuff in 1/m, it only takes 1/mm

%%

% figure();
% hold on;
% %PLT=reshape(ktraj,);
% for i = 1:238
%     %for j = 1:128
%     plot((ktraj(i,:,1)));
%     %old;
%     pause(0.01);
%     %end
%     
% end



%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
fname = ['Rad3D_FID_v1_nospoil_ noTEdel_noTRdel', num2str(radp.Ns),'_',num2str(Nx),'_',num2str(Nz),'_',num2str(rf_dur) '.seq'];
seq.write(fname)
if(wr_traj)
     uisave('ktrajs', 'Ktraj' ); % specify the name and directory to save the k-space trajectory
end
%%
% % Display the first few lines of the output file
% s=fileread('Rad2D_FID.seq');
% disp(s(1:309))
