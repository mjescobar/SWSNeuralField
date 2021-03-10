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
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;

%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=910;
previousTime=6;
samplingFrequency=100;

%Parameters of stimulation
stimFrequency=0.85;
stimAmplitude=1;
stimPulseDuration=0.1;
endTime=905;
startTime=5;
noiseSD=1/sqrt(2)*1e-4;
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
stimMode=0;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);
plasticityFrecuency=1/60;
rateHoldLearn=2;

%Presolver
config=core.Config();
config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
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
stimFrequency=[0.85];
stimPulseDuration=[0.1111,0.0909,0.05,0.011];
seeds=1:5;
interpulsetime=0;
npulses=[3,4];
lenVar1=length(stimFrequency);
lenVar2=length(stimFrequencySD);
lenVar3=length(seeds);
lenVar4=length(stimPulseDuration);
lenVar5=length(npulses);
numSims=lenVar1*lenVar2*lenVar3*lenVar4*lenVar5;
simulationSpace=[lenVar1,lenVar2,lenVar3,lenVar4,lenVar5];

%%%Simulation 1

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);
preseed=1;
simulationTime=tic();
[steadyV,steadyQ,steadyPhi,steadyPreviousV,current,currentQ,currentPhi,delay]=simulator.presolve(stimulator,preseed);
fprintf('Presolver commplete \n')
toc(simulationTime);

parpool(36);
parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3,idxVar4,idxVar5]=ind2sub(simulationSpace,iSim);
	config=core.Config();
	stimAmplitude=sqrt(40/(stimPulseDuration(idxVar4)*npulses(idxVar5)))*sqrt(3);
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
	endTime=ceil(180/stimFrequency(idxVar1)+5); %30 multiple-pulses.
	config=config.setStimParams(stimFrequency(idxVar1),stimAmplitude,stimPulseDuration(idxVar4),startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD(idxVar2),sigmaE,sigmaI,stimX,stimY,shapeNoiseColor,npulses(idxVar5),interpulsetime);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	fprintf('%d\n',stimulator.npulses);
	fprintf('%.4f\n',stimulator.interpulsetime);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	simulationPoint=['SWSlongE-seed-',num2str(seeds(idxVar3)),'-npulses-',num2str(npulses(idxVar5))];
	filename=core.util.buildFilename(simulationPoint,stimShape,stimFrequency(idxVar1), stimAmplitude, stimPulseDuration(idxVar4), stimAmplitudeSD, stimFrequencySD(idxVar2), plasticity_onoff,plasticityFrecuency, rateHoldLearn);
	monitor=core.Monitor(config,['timeseries/multiples/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solve(seeds(idxVar3),steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ)
	fprintf('End simulation of %s \n',filename);
end

fprintf('####THE END####\n');





