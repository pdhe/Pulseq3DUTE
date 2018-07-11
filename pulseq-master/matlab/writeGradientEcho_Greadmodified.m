seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=256; Ny=256;     % Define FOV and resolution

% set system limits
sys = mr.opts('rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6);

% Create 20 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(20*pi/180,'Duration',4e-3,...
    'SliceThickness',5e-3,'apodization',0.5,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',3e-3); %%%%%%%%%%%%%%%%%%%%%%%%
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime); % what is gx.flatTime???
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',2e-3);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',2e-3);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;  

% Calculate timing   % this part it different from the plot. 
delayTE=10e-3 - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2 ...
    - mr.calcDuration(gx)/2;
delayTR=1000e-3 - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
    - mr.calcDuration(gx) - delayTE;

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    seq.addBlock(rf,gz); % what about Gz????? phaseAreas is not a function?
    gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',2e-3); % we have to change gy very time. 
    seq.addBlock(gxPre,gyPre,gzReph);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx,adc);
    seq.addBlock(mr.makeDelay(delayTR))
end

seq.write('gre_Greadmodified_Pei.seq')       % Write to pulseq file
seq.plot()