clear
addpath('./core');
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
%initialStrengths=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];
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
stimAmplitude=12.248*2;
stimPulseDuration=0.1;
endTime=65;
startTime=5;
noiseSD=1/sqrt(2)*1e-4;;
noiseMean=1;
noiseColor=1;
stimShape=13;
targetPhase=0;
stimAmplitudeSD=0;
stimFrequencySD=0;
sigmaE=1;
sigmaI=2;
stimX=7;
stimY=7;
shapeNoiseColor=0;
plasticity_onoff=1;
stimMode=1;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);
plasticityFrecuency=1/60;
rateHoldLearn=2;


config=core.Config();
config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);

stimulator=core.Stimulator(stimMode,config);
model=core.Model(config);
integrator=core.Integrator(config,model);
plasticity=core.Plasticity(config);
monitor=core.Monitor(config,'timeseries/Sleep-closed',{"Phi:E:all"});
simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);

stimFrequency=[0.85,1.0];
targetPhase=(0:45:90)*pi/180;%degrees
seeds=1:5;

lenVar1=length(targetPhase);
lenVar2=length(seeds);
lenVar3=length(stimFrequency);
numSims=lenVar1*lenVar2*lenVar3;
simulationSpace=[lenVar1,lenVar2,lenVar3];



%%%Simulation 1

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);
preseed=1;
simulationTime=tic();
[steadyV,steadyQ,steadyPhi,steadyPreviousV,current,currentQ,currentPhi,delay]=simulator.presolve(stimulator,preseed);
fprintf('Presolver commplete \n')
toc(simulationTime);


%parpool creation
parpool(36);

parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3]=ind2sub(simulationSpace,iSim);
	config=core.Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
	endTime=ceil(30/stimFrequency(idxVar3)+5); %30 pulses.
	config=config.setStimParams(stimFrequency(idxVar3),stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase(idxVar1),stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	simulationPoint=['SWS-seed-',num2str(seeds(idxVar2)),'-phase-',sprintf('%.2f',targetPhase(idxVar1))];
	filename=core.util.buildFilename(simulationPoint,stimShape,stimFrequency(idxVar3), stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff,plasticityFrecuency, rateHoldLearn);
	monitor=core.Monitor(config,['timeseries/phase/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	monitor=monitor.createFilePhase(['timeseries/phase/',filename]);
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solve(seeds(idxVar2),steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ);
	fprintf('End simulation of %s \n',filename);
end

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);
fprintf('#### The END ####');
