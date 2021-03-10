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
%N3
initialStrengths=[6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.58e-3];
%N2%initialStrengths=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];
%initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
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
noiseColor=1;
stimShape=0;
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



%seeds=[234,345,456,567,678,789,901,987,876,765,654,543,432,321];
%letter=['g','h','i','j','k','l','m','n','o','p','q','r','s','t'];
%lenVar2=length(seeds);

%%%%%
%Simulations preparation
%Variable parameters
noiseMean=[0,1];
noiseColor=[1,2];
stimAmplitude=0;
seed=10;
lenVar1=length(noiseMean);
lenVar2=length(noiseColor);
numSims=lenVar1*lenVar2;
simulationSpace=[lenVar1,lenVar2];

%%%Simulation 1

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);



%parpool creation
parpool(4);

parfor iSim=1:numSims
	[idxVar1,idxVar2]=ind2sub(simulationSpace,iSim);
	config=core.Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean(idxVar1),noiseColor(idxVar2),stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	simulationPoint=['N3-noisemean-',num2str(noiseMean(idxVar1)),'-ncolor-',num2str(noiseColor(idxVar2))];
	filename=core.util.buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff,plasticityFrecuency, rateHoldLearn);
	monitor=core.Monitor(config,['timeseries/',filename],{"Phi:E:all","Phi:N:all"});
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	preseed=10;
	[steadyV,steadyQ,steadyPhi,steadyPreviousV,current,currentQ,currentPhi,delay]=simulator.presolve(stimulator,preseed);
	fprintf('Presolver commplete \n')
	simulator=simulator.solve(seed,steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ)
	fprintf('End simulation of %s \n',filename);
end

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);




fprintf('####THE END####\n');


