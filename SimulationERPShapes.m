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
%AAS Model
tau_v=10.0; 
tau_m=10.0;
nu_vc=-2.9e-3;
nu_vh=1e6;
nu_vm=-2.1e-3;
nu_mv=-1.8e-3;
xi=45*3600;
mu=4.4e-9;
c0=4.5;
Cphase=2*pi/3;
Ach=1.3e-3;
AASconnections=1;
%Plasticity %Not used
ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45]*1e-5;
tp=0.01;

%Variables
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
seeds=1;

%Variables behavior
noiseMean=1;
noiseSD=0;
%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=10;
previousTime=10;
samplingFrequency=100;

%delete any previous pool
%poolobj = gcp('nocreate');
%delete(poolobj);
%parpool(8);

stimShapes=[0,1,2,3,4,13];
stimAmplitudes=[10,17.32,17.32,17.32,15.897,12.248]*2*sqrt(2);
stimFrequency=0.2;
stimPulseDuration=0.05;
endTime=10;
startTime=2;
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


numSims=6;
simulationSpace=[6];

for iSim=1:numSims
	[idxVar1]=ind2sub(simulationSpace,iSim);
	%Solver with presolve
	stimAmplitude=stimAmplitudes(idxVar1);
	stimShape=stimShapes(idxVar1);
	%Configuration object
	config=Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
	config=config.setModelAASParams(nu_vm,nu_mv,nu_vc,nu_vh,xi,mu,c0,tau_m,tau_v,Cphase,Ach,AASconnections);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(ratesA,ratesB,tp);

	%Modules
	stimulator=Stimulator(stimMode,config);
	model=Model(config);
	modelAAS=ModelAAS(config);
	integrator=Integrator(config,model,modelAAS);
	plasticity=Plasticity(config);
	filename=sprintf('ShapePulse-%d',idxVar1);
	monitor=Monitor(config,['timeseries/shapes/',filename],{"Phi:E:all","Phi:N:all"});
	simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solveWithPresolve(seeds);
	fprintf('End simulation of %s \n',filename);
end %End main for

%delete any previous pool
%poolobj = gcp('nocreate');
%delete(poolobj);

fprintf('####THE END####\n');
