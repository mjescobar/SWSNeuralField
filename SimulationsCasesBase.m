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
initialChangeRates=[1,-2,3,4,5,6,-7,8]*1e-3; %not used


%Variables
nu1=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];
nu2=[6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.6e-3];
difference_nus=nu2-nu1;
phins=[0,1,2];
stds=[-1,pi/2,pi];
seeds=1:20;

%Variables behavior
samples_nu=5;
samples_phi=2;
samples_std=1;
seed_n=length(seeds);
nus=zeros(samples_nu,8);
for i=1:samples_nu
	nus(i,:)=nu1;
	nus(i,:)=nu1+1/samples_nu*(i-1)*difference_nus;
end
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
parpool(36);

stimFrequency=1;
stimAmplitude=0;
stimPulseDuration=0;
endTime=1;
startTime=0;
noiseColor=2; %White noise: 1 more fast simulation
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

numSims=samples_nu*samples_phi*samples_std*seed_n;
simulationSpace=[samples_nu,samples_phi,samples_std,seed_n];

parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3,idxVar4]=ind2sub(simulationSpace,iSim);
	%Solver with presolve
	initialStrengths=nus(idxVar1,:);
	noiseMean=phins(idxVar2);
	noiseSD=stds(idxVar3+1);
	config=core.Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);
	stimulator=core.Stimulator(stimMode,config);
	model=core.Model(config);
	integrator=core.Integrator(config,model);
	plasticity=core.Plasticity(config);
	filename=sprintf('NU%d-PHIN%d-STD%d-seed%d-pink',idxVar1,idxVar2,idxVar3+1,idxVar4);
	monitor=core.Monitor(config,['timeseries/baselines/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solveWithPresolve(seeds(idxVar4))
	fprintf('End simulation of %s \n',filename);
end %End main for

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);

fprintf('####THE END####\n');
