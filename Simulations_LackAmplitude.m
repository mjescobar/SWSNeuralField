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
stimFrequency=0.2;
stimAmplitude=1;
stimPulseDuration=0.1;
endTime=65;
startTime=5;
noiseSD=-1;
noiseMean=1;
noiseColor=2;
stimShape=2;
targetPhase=0;
stimAmplitudeSD=0;
stimFrequencySD=0.002;
sigmaE=1;
sigmaI=2;
stimX=7;
stimY=7;
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



stimulator=core.Stimulator(stimMode,config);
model=core.Model(config);
integrator=core.Integrator(config,model);
plasticity=core.Plasticity(config);
%filename=core.util.buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff,plasticityFrecuency, rateHoldLearn)
monitor=core.Monitor(config,'presolver',{"Phi:E:all"});
simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);

%seeds=[234,345,456,567,678,789,901,987,876,765,654,543,432,321];
%letter=['g','h','i','j','k','l','m','n','o','p','q','r','s','t'];
%lenVar2=length(seeds);

%%%%%
%Simulations preparation
%Variable parameters
stimShape=2;
stimAmplitudevar=[1];
stimPulseDuration=[0.01];
stimFrequencySD=[0,100];
seeds=1:20;
lenVar1=length(stimAmplitudevar);
lenVar2=length(stimFrequencySD);
lenVar3=length(seeds);
lenVar4=length(stimPulseDuration);
numSims=lenVar1*lenVar2*lenVar3*lenVar4;
simulationSpace=[lenVar1,lenVar2,lenVar3,lenVar4];

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
parpool(4);

parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3,idxVar4]=ind2sub(simulationSpace,iSim);
	stimAmplitude=173.2*stimAmplitudevar(idxVar1);
	config=core.Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration(idxVar4),startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD(idxVar2),sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	simulationPoint=['SWS-seed-',num2str(seeds(idxVar3))];
	filename=core.util.buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration(idxVar4), stimAmplitudeSD, stimFrequencySD(idxVar2), plasticity_onoff,plasticityFrecuency, rateHoldLearn);
	monitor=core.Monitor(config,['timeseries/amplitude/',filename],{"Phi:E:all","Phi:N:all"});
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solve(seeds(idxVar3),steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ)
	fprintf('End simulation of %s \n',filename);
end

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);




fprintf('####THE END####\n');


