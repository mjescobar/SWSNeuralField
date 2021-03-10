clear

%Parameters of the model
%Constants
Qmax=340;
alpha=45;
beta=186;
gamma=116;
range=0.086;
sigma_rho=0.0038;
theta=0.01292;
t0=0.085;
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
initialChangeRates=[1,-2,3,4,5,6,-7,8]*1e-3; %not used



phins=[1];
stds=[1/sqrt(2)*1e-4];
seeds=1:20;

%Variables behavior
samples_phi=1;
samples_std=1;
seed_n=length(seeds);
%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=70;
previousTime=6;
samplingFrequency=100;

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);
parpool(4);

stimFrequency=1;
stimAmplitude=0;
stimPulseDuration=0;
endTime=1;
startTime=0;
noiseColor=1; %White noise: 1 more fast simulation
stimShape=0;
targetPhase=0;
stimAmplitudeSD=0;
stimFrequencySD=0;
sigmaE=1;
sigmaI=2;
stimX=7;
stimY=7;
shapeNoiseColor=0;
plasticity_onoff=0; %Activate with 2
stimMode=0;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);
plasticityFrecuency=1/60; %not used
rateHoldLearn=2;% not used

numSims=samples_phi*samples_std*seed_n;
simulationSpace=[samples_phi,samples_std,seed_n];

parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3]=ind2sub(simulationSpace,iSim);
	%Solver with presolve
	noiseMean=phins(idxVar1);
	noiseSD=stds(idxVar2);
	config=core.Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	filename=sprintf('SWS-BaselineWhiteNoise-seed%d',idxVar3);
	monitor=core.Monitor(config,['timeseries/baselines/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solveWithPresolve(seeds(idxVar3))
	fprintf('End simulation of %s \n',filename);
end %End main for

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);

fprintf('####THE END####\n');
