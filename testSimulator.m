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
%initialStrengths=[3.06e-3, -3.24e-3, 0.92e-3, 0.26e-3, 2.88e-3, 4.73e-3, -1.95e-3, 2.70e-3];
initialStrengths=[5.54 -5.65 1.53 0.286 1.12 2.67 -1.73 9.22]*1e-3;
initialChangeRates=[1,-2,3,4,5,6,-7,8]*1e-3;

%Parameters of the simulation
Nx=16;
Ny=16;
Lx=0.5;
Ly=0.5;
h=1e-4;
finalTime=5400;
previousTime=5;
samplingFrequency=100;

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
rateHoldLearn=2;


config=core.Config();
config=config.setTimeParams(h,finalTime,previousTime,samplingFrequency);
config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,initialChangeRates,plasticity_onoff);
config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor);
config=config.setPlasticityParams(plasticityFrecuency,rateHoldLearn);

stimulator=core.Stimulator(stimMode,config);
model=core.Model(config);
integrator=core.Integrator(config,model);
plasticity=core.Plasticity(config);
monitor=core.Monitor(config,'timeseries/Sleep-90min',{"Phi:E:all","Phi:N:all"});
simulator=core.Simulator(config,model,integrator,stimulator,monitor,plasticity);

seed=1;
simulationTime=tic();
simulator=simulator.solveWithPresolve(seed);
toc(simulationTime);

% tic()
% for n=1:10000
% [backgroundNoise,noiseFlags]=stimulator.generate_noise(noiseFlags);
% [stim(:,n),marker,stimFlags]=stimulator.generate_stim(n,stimFlags,backgroundNoise);
% end
% toc()
