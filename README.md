# SWSNeuralField
Matlab implementation of the **Neural Field Theory** model.

The intention of this implementation is to build a testbed of stimulation protocols on this large-scale brain model allowing also closed-loop stimulation.


Model theory and slow-wave sleep tunning details are in:
>Selection of stimulus parameters for enhancing slow wave sleep events with a Neural-field theory thalamocortical computational model
>(https://www.biorxiv.org/content/10.1101/2021.02.03.429528v1.article-metrics)

More details of the model could be found here:

> Sanz-Leon, P., Robinson, PA., Knock, SA., Drysdale, PM., Abeysuriya, RG., Fung, PK., Rennie, CJ and Zhao, XL.
>NFTsim: Theory and Simulation of Multiscale Neural Field Dynamics
>PLoS Computational Biology (2018).

And a C++ implementation of the main model is in the repository:

[https://github.com/BrainDynamicsUSYD/nftsim](https://github.com/BrainDynamicsUSYD/nftsim)

## Installation

Clone or download this repository.

``` git clone https://github.com/mjescobar/SWSNeuralField.git ```

To avoid the data collected for *'Selection of stimulus parameters for enhancing slow wavesleep events with a Neural-field theory thalamocortical model'* clone with the command __(requires git >2.30.0)__:

``` git clone --depth 1 --sparse --filter=blob:none https://github.com/mjescobar/SWSNeuralField.git ```

``` cd SWSNeuralField ```

``` git sparse-checkout init --cone ```

``` git spatse-checkout set core ```

### Requirements

Simulations tested on Matlab R2017b on Ubuntu 18.04 (Server and Desktop) and R2020b on Ubuntu 20.04.

Analysis code tested on Python >3.6.9. Package requirements in 'analyze/requirements.txt' 

The requirements could be installed from inside that folder with (preferably also inside a python virtual environment):

```pip install -r requirements.txt```

or (if python2 is also installed in the system): 

``` pip3 install -r requirements.txt ```

## Use
### Simulation
Create an folder for output files (example: _'./timeseries'_) in the SWSNeuralFiled main directory

Run 

``` matlab -nodisplay -r "run('./<StimuationFile.m>');exit;"& ```

### Configure simulation
The file *testSimulation.m* shows a demo configuration of simulation configuration.

In brief, there is the need to set the model parameters. 

Then, instatiate the configuration object and pass to it the parameters:

``` config=Config(); ```

- Time parameters

``` config.setTimeParams(timeSimulationStep,finalTime,previousStorageTime,samplingFrequency);  ```

- Neural Field model parameters

``` config=config.setModelParams(Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff); ```

- Stimulation parameters

``` config=config.setStimParams(stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor); ```

- Plasticity parameters (_in development, but already required now_)

``` config=config.setPlasticityParams(ratesA,ratesB,tp); ```

After, the simulation modules need to be instatiate:

- Stimulator

``` stimulator=Stimulator(stimMode,config); ```

_stimMode_ could be 0: Open loop , 1: Closed loop.

- Model

``` model=Model(config); ```

- Integrator

``` integrator=Integrator(config,model); ```

- Plasticity (_in development, but already required now_)

``` plasticity = Plasticity(config); ```

- Monitor

``` monitor=Monitor(config,['timeseries/',filename],{"Phi:E:all","Phi:N:all"}); ```

(If _stimMode_==1, you can also include ``` monitor=monitor.createFilePhase(['timeseries',filename]);```)


The outputs syntaxis is: 
``` "_Measurement_:_Population_:_indexes_" ```  inside a cell's block of Matlab "_{}_".

Where **_Measurement_** could be **Phi**: propagation field, **V**: soma voltage, **Q**: firing response.

**_Population_** is one of the populations, **E**: excitatory cortex, **I**: inhibitory cortex, **R**: reticular nucleus, **S**: relay nuclei, **N**: input.

**_indexes_** could be **'all'** to include all the nodes (Ny x Nx) or a list of node indexes inside brackets [].

Finally **_filename_** is the name of output files before the suffix. There is an aditional method to crate filenames from stimulation parameters:

``` filename=util.buildFilename(prefix,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff);```

where _prefix_ is any string, preferably indicating something about the model or time parameters (example: 'SWS-90min')

There could be two or three output files (if _stimMode_==1). The files' format is plain text with number format '%.9e', and can be read by any .CSV application or library.

-_filename_-Phi.txt: timeseries of *n x 4* nodes (column division: tabular).

-_filename_-Stim.txt: timeseries of *n* input nodes, time marker and stimulation markers (column division: tabular).

-_filename_-Phase.txt: online filtered signal, online envelope, detected phase (column division: comma).

- Simulator

``` simulator=Simulator(config,model,integrator,stimulator,monitor,plasticity); ```

The simulation could be performed by three functions: *presolve*, *solve*, and *solveWithPresolve*.

*presolve* performs the simulation for the indicated _previousStorageTime_ to obtain initial conditions different than zero. Note that the _stimulator_ is requeired to get the mean and standard deviation of the background noise and return encapuslated the mean value of x(r,t).

``` [steadyV,steadyQ,steadyPhi,steadyPreviousV,current,currentQ,currentPhi,delay]=simulator.presolve(stimulator,preseed); ```

*solve* performs the simulation storing the results at the indicated _samplingFrequency_.

``` simulator.solve(seed,steadyV,steadyQ,steadyPhi,steadyPreviousV,currentQ) ```
	
*solveWithPresolve* performs the entire simulation. Not recommended for multiple simulations by the overcharge of the presolver, with not stored-data.

``` simulator.solveWithPresolve(seed) ```

### Analysis

#### Spectrum analysis

Load simulation data functions and spectrum analysis are in _'spectrum'_ directory.

#### Detection of events

Detection of SWS events functions are located in _'eventsDetection'_ directory.

#### Plot

Auxiliar functions for plots are in _'plotting'_ subdirectory

#### Collected Data

Data from SWS stimulation with different pulse duration, pulse energy, frequencies, and target-phase of closed-loop stimulation are stored in _'collectedData'_ (313 MB).
