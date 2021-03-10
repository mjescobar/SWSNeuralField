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
Cphase=2*pi/3;
Ach=1.3e-3;
AASconnections=1;
%Variables
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
previousTime=10;
samplingFrequency=100;

%Parameters of stimulation
stimFrequency=0.001;
stimAmplitude=0;
stimPulseDuration=0.1;
endTime=6;
startTime=5;
noiseSD=1/sqrt(2)*1e-4;
noiseMean=0;
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
plasticity_onoff=2; %2:ON
stimMode=0;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);


%Presolver
config=Config();
config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
config=config.setModelAASParams(nu_vm,nu_mv,nu_vc,nu_vh,xi,mu,c0,tau_m,tau_v,Cphase,Ach,AASconnections);
config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
config=config.setPlasticityParams(ratesA,ratesB,tp);



stimulator=Stimulator(stimMode,config);
model=Model(config);
modelAAS=ModelAAS(config);
integrator=Integrator(config,model,modelAAS);
plasticity=Plasticity(config);
plasticity=plasticity.init();
monitor=Monitor(config,'presolver',{"Phi:E:all"});
simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);

%seeds=[234,345,456,567,678,789,901,987,876,765,654,543,432,321];
%letter=['g','h','i','j','k','l','m','n','o','p','q','r','s','t'];
%lenVar2=length(seeds);

%%%%%
%Simulations preparation
%Variable parameters
stimShape=0;
stimPulseDuration=0.01;
stimFrequencySD=0;
seeds=1;
AASconnections=1;
lenVar1=length(AASconnections);
lenVar2=length(seeds);
numSims=lenVar1*lenVar2;
simulationSpace=[lenVar1,lenVar2];
%%%Simulation 1

%delete any previous pool
%poolobj = gcp('nocreate');
%delete(poolobj);
preseed=1;
simulationTime=tic();
[steadyV,steadyQ,steadyPhi,steadyPreviousV,Vaas,Qaas,current,currentQ,currentPhi,delay]=simulator.presolveAAS(stimulator,preseed);
fprintf('Presolver commplete \n')
toc(simulationTime);

%parpool creation
%parpool(6);

for iSim=1:numSims
	[idxVar1,idxVar2]=ind2sub(simulationSpace,iSim);
	config=Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
	config=config.setModelAASParams(nu_vm,nu_mv,nu_vc,nu_vh,xi,mu,c0,tau_m,tau_v,Cphase,Ach,AASconnections);
	config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(ratesA,ratesB,tp);
	stimulator=Stimulator(stimMode,config);
	model=Model(config);
	modelAAS=ModelAAS(config);
	integrator=Integrator(config,model,modelAAS);
	plasticity=Plasticity(config);
	plasticity=plasticity.init();
	simulationPoint=['SWSAASPhi0',num2str(AASconnections),'-seed-',num2str(seeds(idxVar2))];
	filename=util.buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff);
	monitor=Monitor(config,['timeseries/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	monitor=monitor.createFilePhase(['timeseries/',filename]);
	monitor=monitor.createFileStrengths(['timeseries/',filename]);
	simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solveAAS(seeds(idxVar2),steadyV,steadyQ,steadyPhi,steadyPreviousV,Vaas,Qaas,currentQ)
	fprintf('End simulation of %s \n',filename);
end

%delete any previous pool
%poolobj = gcp('nocreate');
%delete(poolobj);




fprintf('####THE END####\n');


