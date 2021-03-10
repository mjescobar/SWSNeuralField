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
ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45]*1e-5;
tp=0.01;
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
stimAmplitude=17.32*2;
%stimPulseDuration=0.10; %DONE

stimPulseDuration=0.05;

endTime=905;
startTime=5;
noiseSD=1/sqrt(2)*1e-4;;
noiseMean=1;
noiseColor=1;
stimShape=2;
targetPhase=0;
stimAmplitudeSD=0;
stimFrequencySD=0;
sigmaE=1;
sigmaI=2;
stimX=7;
stimY=7;
shapeNoiseColor=0;
plasticity_onoff=0;
stimMode=1;
stimFlags=zeros(1,10);
noiseFlags=zeros(1,3);



config=Config();
config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
config=config.setPlasticityParams(ratesA,ratesB,tp);

stimulator=Stimulator(stimMode,config);
model=Model(config);
modelAAS=ModelAAS(config);
integrator=Integrator(config,model,modelAAS);
plasticity=Plasticity(config);
monitor=Monitor(config,'timeseries/Sleep-closed',{"Phi:E:all"});
simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);

stimFrequency=[0.2,0.5,0.85];
stimAmplitude=[17.31*2*sqrt(2)];
stimPulseDuration=[0.05];
targetPhase=(0:45:90)*pi/180;%degrees
seeds=1:5;

lenVar1=length(targetPhase);
lenVar2=length(seeds);
lenVar3=length(stimFrequency);
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
parpool(36);

parfor iSim=1:numSims
	[idxVar1,idxVar2,idxVar3,idxVar4]=ind2sub(simulationSpace,iSim);
	config=Config();
	config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
	config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff);
	endTime=ceil(180/stimFrequency(idxVar3)+5); %30 pulses.
	config=config.setStimParams(stimFrequency(idxVar3),stimAmplitude(idxVar4),stimPulseDuration(idxVar4),startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase(idxVar1),stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
	config=config.setPlasticityParams(ratesA,ratesB,tp);
	stimulator=Stimulator(stimMode,config);
	model=Model(config);
	modelAAS=ModelAAS(config);
	integrator=Integrator(config,model,modelAAS);
	plasticity=Plasticity(config);
	simulationPoint=['SWSlong-seed-',num2str(seeds(idxVar2)),'-phase-',sprintf('%.2f',targetPhase(idxVar1))];
	filename=util.buildFilename(simulationPoint,stimShape,stimFrequency(idxVar3), stimAmplitude(idxVar4), stimPulseDuration(idxVar4), stimAmplitudeSD, stimFrequencySD, plasticity_onoff);
	monitor=Monitor(config,['timeseries/phase/',filename],{"Phi:E:all","Phi:N:all","Phi:R:135","Phi:S:135"});
	monitor=monitor.createFilePhase(['timeseries/phase/',filename]);
	simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);
	fprintf('Start simulation of %s \n',filename);
	simulator=simulator.solve(seeds(idxVar2),steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ);
	fprintf('End simulation of %s \n',filename);
end

%delete any previous pool
poolobj = gcp('nocreate');
delete(poolobj);
fprintf('#### The END ####');
