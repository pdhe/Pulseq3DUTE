seq=mr.Sequence();              % Create a new sequence object
Nx=256;                         % number of voxels of a line in k space?
Nrep=1;
addpath(genpath('.'));          % why do we need to addpath? 
%% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',10e-3);

%% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',3.2e-3);
delayTE=20e-3;
delayTR=100e-3;
                                % where is our space encoding process?
%% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf);                   %rf
    seq.addBlock(mr.makeDelay(delayTE));%wait
    seq.addBlock(adc);                  %receive
    seq.addBlock(mr.makeDelay(delayTR)) %wait
end
%%
seq.write('fid.seq')       % Write to pulseq file
seq.plot('TimeRange',[0 100e-3]);