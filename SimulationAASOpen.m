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
tau_v=10.0; 
tau_m=10.0;
nu_vc=-2.9e-3;
nu_vh=1e6;
nu_vm=-2.1e-3;
nu_mv=-1.8e-3;
xi=45*3600;
mu=4.4e-9;
c0=4.5;
%Cphase=2*pi/3;
Ach=1.3e-3;
AASconnections=1;
%Variables of connections strengths
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45]*1e-5;
tp=0.01;
%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=1800;
previousTime=60; %1 minute
samplingFrequency=100;

%Parameters of stimulation
stimFrequency=0.85;
stimPulseDuration=0.01;
stimShape=2;
stimAmplitude=17.32*2;
endTime=1500;
startTime=5;
noiseSD=1/sqrt(2)*1e-4;
noiseMean=0;
noiseColor=1;
stimShape=0;
targetPhase=0;
stimAmplitudeSD=0;
sigmaE=1;
sigmaI=2;
stimX=7;
stimY=7;
shapeNoiseColor=0;
plasticity_onoff=2; %2:ON
stimMode=0;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);
AASconnections=1;




%%%%%
%Simulations preparation
%Variable parameters

Cphase=[pi];
seeds=4:5;
stimFrequencySD=100;
lenVar1=length(Cphase);
lenVar2=length(seeds);
lenVar3=length(stimFrequencySD);
numSims=lenVar1*lenVar2*lenVar3;
simulationSpace=[lenVar1,lenVar2,lenVar3];
%%%Simulation 1

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);



%parpool creation
parpool(8);

parfor iSim=1:numSims
	[idxVar1,idxVar2,indxVar3]=ind2sub(simulationSpace,iSim);
	config=Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
	config=config.setModelAASParams(nu_vm,nu_mv,nu_vc,nu_vh,xi,mu,c0,tau_m,tau_v,Cphase(idxVar1),Ach,AASconnections);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD(indxVar3),sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(ratesA,ratesB,tp);
	stimulator=Stimulator(stimMode,config);
	model=Model(config);
	modelAAS=ModelAAS(config);
	integrator=Integrator(config,model,modelAAS);
	plasticity=Plasticity(config);
	plasticity=plasticity.init();
	simulationPoint=['SWSAASPhin0-seed-',num2str(seeds(idxVar2)),'-Cphase',num2str(Cphase(idxVar1))];
	filename=util.buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD(indxVar3), plasticity_onoff);
	monitorpre=Monitor(config,['presolver-',filename],{"Phi:E:1"});
	simulatorpre=Simulator(config,model,modelAAS,integrator,stimulator,monitorpre,plasticity);
	[steadyV,steadyQ,steadyPhi,steadyPreviousV,Vaas,Qaas,current,currentQ,currentPhi,delay]=simulatorpre.presolveAAS(stimulator,seeds(idxVar2));
	fprintf('Completed Presolver of %s \n',filename)
	monitor=Monitor(config,['timeseries/aas/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	monitor=monitor.createFilePhase(['timeseries/aas/',filename]);
	monitor=monitor.createFileStrengths(['timeseries/aas/',filename]);
	simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solveAAS(seeds(idxVar2),steadyV,steadyQ,steadyPhi,steadyPreviousV,Vaas,Qaas,currentQ)
	fprintf('End simulation of %s \n',filename);
end

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);




fprintf('####THE END####\n');


