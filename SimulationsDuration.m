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
%Variables
%N2
%initialStrengths=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];
%66% N2 to N3 (Abysuriya,2015)
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
initialChangeRates=[1,-2,3,4,5,6,-7,8]*1e-3;

%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=70;
previousTime=6;
samplingFrequency=100;

%Parameters of stimulation
stimFrequency=0.85;
stimAmplitude=0;
stimPulseDuration=0.05;
endTime=65;
startTime=5;
%BackgroundNoise
noiseSD=-1;
noiseMean=1;
noiseColor=1;
stimShape=0;
targetPhase=0;
sigmaE=2;
sigmaI=3;
stimX=11;
stimY=13;
%Stim noise
stimAmplitudeSD=0;
stimFrequencySD=0;
shapeNoiseColor=0;
plasticity_onoff=1;
stimMode=0;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);
plasticityFrecuency=1/60;
rateHoldLearn=2;

%Presolver
config=core.Config();
config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);

simulationPoint='SWS-mean1';

stimulator=core.Stimulator(stimMode,config);
model=core.Model(config);
integrator=core.Integrator(config,model);
plasticity=core.Plasticity(config);
%filename=core.util.buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff,plasticityFrecuency, rateHoldLearn)
monitor=core.Monitor(config,'presolver',{"Phi:E:all"});
simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);

seed=1;
simulationTime=tic();
[steadyV,steadyQ,steadyPhi,steadyPreviousV,current,currentQ,currentPhi,delay]=simulator.presolve(stimulator,seed);
fprintf('Presolver commplete \n')
toc(simulationTime);

%%%%%
%Simulations preparation
%Variable parameters
stimFrequency=[0.5, 0.6, 0.7, 0.8, 0.9, 1];
stimAmplitude=[1, 5, 10, 20, 30, 40, 50];
stimPulseDuration=[0.05, 0.1, 0.15, 0.2, 0.25, 0.3];
stimFrequencySD=[0, 0.01];
stimShape=[0,1,2,3,4,5,6,7,8,9,10,11];

lenVar1=length(stimFrequency);
lenVar2=length(stimAmplitude);
lenVar3=length(stimPulseDuration);
lenVar4=length(stimShape);
lenVar5=length(stimFrequencySD);

numSims=lenVar1*lenVar2*lenVar3*lenVar4*lenVar5;
simulationSpace=[lenVar1,lenVar2,lenVar3,lenVar4,lenVar5];

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);

%parpool creation
parpool(30);
parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3,idxVar4,idxVar5]=ind2sub(simulationSpace,iSim);
	config=core.Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
	config=config.setStimParams(stimFrequency(idxVar1),stimAmplitude(idxVar2),stimPulseDuration(idxVar3),startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape(idxVar4),targetPhase,stimAmplitudeSD,stimFrequencySD(idxVar5)*stimFrequency(idxVar1),sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	filename=core.util.buildFilename(simulationPoint,stimShape(idxVar4),stimFrequency(idxVar1), stimAmplitude(idxVar2), stimPulseDuration(idxVar3), stimAmplitudeSD, stimFrequencySD(idxVar5)*stimFrequency(idxVar1), plasticity_onoff,plasticityFrecuency, rateHoldLearn);
	monitor=core.Monitor(config,['timeseries/',filename],{"Phi:E:all"});
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solve(steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ)
	fprintf('End simulation of %s \n',filename);
end
fprintf('####THE END####\n');


