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
%AAS Model %Not used
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
%initialStrengths=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
%Plasticity
ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45]*1e-5;
tp=0.01;

%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=5400; %seconds
previousTime=5; %seconds
samplingFrequency=100; %secondds

%Parameters of stimulation
stimFrequency=0.8;
stimAmplitude=2;
stimPulseDuration=0.05;
endTime=60;
startTime=30;
noiseSD=-1;
noiseMean=0;
noiseColor=1;
stimShape=0;
targetPhase=0;
stimAmplitudeSD=0;
stimFrequencySD=100;
sigmaE=1;
sigmaI=2;
stimX=7;
stimY=7;
shapeNoiseColor=0;
plasticity_onoff=0;
stimMode=0;
stimFlags=zeros(1,6);
noiseFlags=zeros(1,3);
plasticityFrecuency=1/60;

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
monitor=Monitor(config,'timeseries/Sleep-90min',{"Phi:E:all","Phi:N:all"});
simulator=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity);

%Run simulation
seed=1;
simulationTime=tic();
simulator=simulator.solveWithPresolve(seed);
toc(simulationTime);

%Test stimulator
% tic()
% for n=1:10000
% [backgroundNoise,noiseFlags]=stimulator.generate_noise(noiseFlags);
% [stim(:,n),marker,stimFlags]=stimulator.generate_stim(n,stimFlags,backgroundNoise);
% end
% toc()
