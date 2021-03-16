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
nu1=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];%Spindles
nu2=[6.81e-3, -6.89e-3, 1.85e-3, 0.3e-3, 0.21e-3, 1.61e-3, -1.62e-3, 12.6e-3];%N3
difference_nus=nu2-nu1;
seeds=1;

%Variables behavior
samples_nu=5;
phin=[0,0,0,1,2];
stdn=1/sqrt(2)*1e-4;
seed_n=length(seeds);
nus=zeros(samples_nu,8);
nus(1,:)=nu1;
nus(2,:)=nu2;
nus(3,:)=nu1+2/3*difference_nus;
nus(4,:)=nu1+2/3*difference_nus;
nus(5,:)=nu1+2/3*difference_nus;
%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=20;
previousTime=10;
samplingFrequency=100;

%delete any previous pool
poolobj = gcp('nocreate');
%delete(poolobj);
%parpool(8);

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


numSims=samples_nu;
simulationSpace=[samples_nu];

for iSim=1:numSims
	[idxVar1]=ind2sub(simulationSpace,iSim);
	%Solver with presolve
	initialStrengths=nus(idxVar1,:);
	noiseMean=phin(idxVar1);
	noiseSD=stdn;
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
	filename=sprintf('Sleep-caso-%d',idxVar1);
	monitor=Monitor(config,['timeseries/baselines/',filename],{"Phi:E:all","Phi:N:all"});
	simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solveWithPresolve(seeds);
	fprintf('End simulation of %s \n',filename);
end %End main for

%delete any previous pool
%poolobj = gcp('nocreate');
%delete(poolobj);

fprintf('####THE END####\n');
