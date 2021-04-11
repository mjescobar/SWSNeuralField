%Config: A matlab to set/store the principal parameters of simulation
classdef Config
    
    properties
        %Time
        h
        finalTime
        previousTime
        samplingFrequency
        samplingLeap
        storageLeap
        storageRows
        %Model
        Nx
        Ny
        Lx
        Ly
        dx
        plasticity_onoff
        alpha
        beta
        Qmax
        sigma_rho
        theta
        gamma
        range
        t0
        initialStrengths
        initialChangeRates
        connectionMatrix
        connectionMatrixDelay
        connectionMatrixNodelay
        plasticityMatrix
        M
        %Stimulation
        stimFrequency
        stimAmplitude
        stimPulseDuration
        stimPulseDurationSamples
        endTime
        startTime
        noiseSD
        noiseMean
        noiseColor
        stimShape
        targetPhase
        stimAmplitudeSD
        stimFrequencySD
        sigmaE
        sigmaI
        stimX
        stimY
        shapeNoiseColor
        stimTimeMarkers
        %Plasticity changes
        plasticityFrequency
		plasticityFlag
        omega
        ratesA
        ratesB
        tp
        %Multiple pulses
        npulses
        interpulsetime
        
        %AAS model
        nu_vm
		nu_mv
		nu_vc
		nu_vh
		xi
		c0
		mu
		tau_m
		tau_v
		Cphase
		Ach   
		AASconnections
		

        
    end
    
    methods
		%Construnctor
        function obj=Config()
            obj.M=5;
            %defaut AAS parameters
            obj.tau_m=10;
            obj.tau_v=10;
            obj.xi=45*3600;
            obj.c0=4.5;
            obj.mu=4.4e-9;
            obj.Ach=1.3e-3;
        end
        %Set time parameters of the simulation (-1 arguments almost always give the default value)
        function obj=setTimeParams(obj,h,finalTime,previousTime,samplingFrequency)
            if h>0
                obj.h=h;
            else
                fprintf('Config: Solver step time must be positive non-zero\n')
                fprintf('Config: Set h=1e-4\n')
                obj.h=1e4;
            end
            if mod(samplingFrequency,1)==0
                dt_save=1/samplingFrequency;
            else
                fprintf('Config: Approximating sampling frequency to the near superior integer\n')
                samplingFrequency=ceil(samplingFrequency);
                dt_save=1/ceil(samplingFrequency);
            end
            if dt_save>h
                obj.samplingFrequency=samplingFrequency;
                obj.samplingLeap=ceil(1/(obj.h*samplingFrequency));
                
            else
                fprintf('Config: Save time must be greather than solver step \n');
                fprintf('Config: Setting sampling frequency to 0.01/h %f: \n',0.01/obj.h);
                obj.samplingFrequency=ceil(0.01/obj.h);
                obj.samplingLeap=ceil(1/(obj.h*samplingFrequency));
            end
            %%% Plasticity and Storage %%%
            %For analytical and plasticity calculations
            %Fixed values
            obj.storageLeap=2*obj.samplingFrequency*obj.samplingLeap;
            obj.storageRows=floor(2*obj.samplingFrequency);
            obj.omega=pi/obj.samplingFrequency:pi/obj.samplingFrequency:pi;
            %%%
            
            %Offset time to save results must be zero or greater than one
            if previousTime>=0
                obj.previousTime=previousTime;
            else
                obj.previousTime=3;
                fprintf('Config: Start time offset must be positive, set =3\n');
            end
            if finalTime>0 && finalTime>obj.previousTime
                obj.finalTime=finalTime;
            else
                fprintf('Config: Final time must be positive non-zero \n');
                fprintf('Config: Final time set to result_start_time+1: %f \n',obj.previousTime+1);
                obj.finalTime=obj.previousTime+1;
            end

        end
            
       %Set time-space parameters of the model (-1 arguments almost always give the default value)
        function obj=setModelParams(obj,Nx,Ny,Lx,Ly,alpha,beta,Qmax,sigma_rho,theta,gamma,range,t0,initialStrengths,plasticity_onoff)

            if Nx>=1 && Nx<50 %At least one node per row and no more than 50 nodes per row
                obj.Nx=Nx;
                fprintf('Config: Set nodes per column (Number of rows) to %d \n',Nx);
            else
                obj.Nx=1;
                fprintf('Config: Set nodes per column to 1 \n');
            end
            if Ny>=1 && Ny<50 %At least one node per column and no more than 50 nodes per column
                obj.Ny=Ny;
                fprintf('Config: Set nodes per row (Number of columns) to %d \n',Nx);
            else
                obj.Ny=1;
                fprintf('Config: Set nodes per row to 1 \n');
            end
            if Lx>0
                obj.Lx=Lx;
            else
                obj.Lx=0.5;
                fprintf('Config: Horizontal length to 0.5 \n');
            end
            if Ly>0
                obj.Ly=Ly;
            else
                obj.Ly=0.5;
                fprintf('Config: Vertical length to 0.5 \n');
            end
            dxp=Lx/Nx;

            if dxp>0
                obj.dx=dxp;
            else
                fprintf('Config: Delta space must be positive non-zero \n');
                fprintf('Config: Set dx=0.01\n');
                obj.dx=0.01;
            end
            
            if alpha<beta
                obj.alpha=alpha;
                obj.beta=beta;
            else
                fprintf('Config: Beta must be grater than alpha \n')
                fprintf('Config: Set default values of sleep alpha=45, beta=186 \n')
                obj.alpha=45;
                obj.beta=186;
            end

            if Qmax>=1
                obj.Qmax=Qmax;
            else
                fprintf('Config: Qmax must be grater or equal than 1 \n')
                fprintf('Config: Set default value Qmax=340 \n')
                obj.Qmax=340;
            end
            if sigma_rho>0
                obj.sigma_rho=sigma_rho;
            else
                fprintf('Config: sigma_rho must be grater than zero \n')
                fprintf('Config: Set default value sigma_rho=0.0038 \n')
                obj.sigma_rho=0.0038;
            end
            if t0>=0
				obj.t0=t0;
			else
				fprintf('Config: t0 can not be negative, set t0=0.0085')
				obj.t0=0.085;
			end
			
            obj.theta=theta;
            obj.gamma=gamma;
            obj.range=range;

            size_ic=size(initialStrengths);
            if size_ic(1)==8 && size_ic(2)==1
                obj.initialStrengths=initialStrengths;
            elseif size_ic(1)==1 && size_ic(2)==8
                obj.initialStrengths=initialStrengths';
            else
                fprintf('Config: A 8x1 or 1x8 vector expected \n')
                fprintf('Config: Set with random vaules in range [-1e-3,1e-3] and with signs constraints \n')
                obj.initialStrengths=rand(8,1)*1e-3;
                obj.initialStrengths(2)=-obj.initialStrengths(2);
                obj.initialStrengths(7)=-obj.initialStrengths(7);
            end
            eirs_connections=util.eirsConnections(obj.initialStrengths);
            [obj.connectionMatrix,obj.connectionMatrixDelay,obj.connectionMatrixNodelay]=util.expandConnections(obj.Nx,obj.Ny,eirs_connections);
            obj.plasticity_onoff=plasticity_onoff;
        end
            
        %Set paaramters of the AAS model
        function obj=setModelAASParams(obj,nu_vm,nu_mv,nu_vc,nu_vh,xi,mu,c0,tau_m,tau_v,Cphase,Ach,AASconnections)
			obj.nu_vm=nu_vm;
			obj.nu_mv=nu_mv;
			obj.nu_vc=nu_vc;
			obj.nu_vh=nu_vh;
			obj.xi=xi;
			obj.c0=c0;
			obj.mu=mu;
			obj.tau_m=tau_m;
			obj.tau_v=tau_v;
			obj.Cphase=Cphase;
			obj.Ach=Ach;
			obj.AASconnections=AASconnections;
		end
        %Set parameters of the stimulator (-1 arguments almost always give the default value)
        function obj=setStimParams(obj,stimFrequency,stimAmplitude,stimPulseDuration,startTime,endTime,noiseSD,noiseMean,noiseColor,stimShape,targetPhase,stimAmplitudeSD,stimFrequencySD,sigmaE,sigmaI,stimX,stimY,shapeNoiseColor,npulses,interpulsetime)
            %Number of multiple pulses
            if nargin>18
				fprintf('Configuring multiple pulses \n')
				fprintf('input npulses: %d\n',npulses)
				fprintf('input interpulsetime: %.3f\n',interpulsetime)
				obj.npulses=npulses;
				obj.interpulsetime=interpulsetime;
				fprintf('Number of pulses: %d\n',obj.npulses)
				fprintf('Interpulse time: %.3f\n',obj.interpulsetime)
			else
				fprintf('Single pulse \n')
				obj.npulses=1;
				obj.interpulsetime=0;
			end
			
            
            if stimFrequency<0.5*obj.samplingFrequency && stimFrequency>0
                obj.stimFrequency=stimFrequency;
            elseif stimFrequency==-1
				%fixed number of pulses
				%For use with phase stimulation 
				obj.stimFrequency=-1;
            else
                fprintf('Config: Stimulation frequency must be lesser than 1/2*sampling frequency \n')
                fprintf('Config: Set with default value stimFrequency=0.5Hz \n')
                obj.stimFrequency=0.5;
            end
            if stimAmplitude>=0
                obj.stimAmplitude=stimAmplitude;
            else
                fprintf('Config: Amplitude must be defined positive \n');
                fprintf('Config: Set with default value stimAmplitude=0 \n');
                obj.stimAmplitude=0;
            end
            if stimPulseDuration>=0
                obj.stimPulseDuration=stimPulseDuration;
                obj.stimPulseDurationSamples=ceil(1/obj.h*obj.stimPulseDuration);
            else
                fprintf('Config: Duration must be defined positive \n');
                fprintf('Config: Set with default value stimPulseDuration=0.01 \n');
                obj.stimPulseDuration=0.01;
                obj.stimPulseDurationSamples=ceil(1/obj.h*obj.stimPulseDuration);
            end
            if startTime>=0 && startTime<=obj.finalTime
                obj.startTime=startTime;
            else
                fprintf('Config: Stimulation start time must be lesser than stimulation total time \n');
                fprintf('Config: Set for all time \n');
                obj.startTime=0;
            end

            if 	endTime>=obj.startTime	&& endTime<=obj.finalTime
                obj.endTime=endTime;
            else
                fprintf('Config: Stimulation end time must be lesser than stimulation total time \n');
                fprintf('Config: Set for all time \n');
                obj.endTime=obj.finalTime;
            end
            if noiseSD>=0
				asd=noiseSD;
				if obj.Nx*obj.Ny>1
                    std_noise=sqrt(2*4*pi^3.*asd.^2/obj.h/obj.dx^2);
                else
                    std_noise=sqrt(2.0*pi.*asd.^2/obj.h);
                end
                obj.noiseSD=std_noise;
            else
                fprintf('Config: Set deafult level of noiseSD from  Abeysuriya2014\n');
                asd=7.071067811865475e-05; %Abeysuriya2014

                if obj.Nx*obj.Ny>1
                    std_noise=sqrt(2*4*pi^3.*asd.^2/obj.h/obj.dx^2);
                else
                    std_noise=sqrt(2.0*pi.*asd.^2/obj.h);
                end
                obj.noiseSD=std_noise;
            end
            if noiseMean>=0
                obj.noiseMean=noiseMean;
            else
                obj.noiseMean=0;
            end
            if noiseColor>=0 && noiseColor<=3
                obj.noiseColor=noiseColor;
                %0 None
                %1 White
                %2 Pink
            else
                fprintf('Config: Set noise color as: White\n');
                obj.noiseColor=1;
            end
            if stimShape>=0 && stimShape<=14
                obj.stimShape=stimShape;
            else
                fprintf('Config: Set rectangular shape\n');
                obj.stimShape=0;
            end

            if targetPhase>=0 && targetPhase<=360
                obj.targetPhase=targetPhase;
            else
                obj.targetPhasetargetPhase=0;
            end
            if stimAmplitudeSD>=0
                obj.stimAmplitudeSD=stimAmplitudeSD;
            else
                fprintf('Config: Deviation of the main amplitude must be defined positive \n');
                obj.stimAmplitudeSD=0;
            end
			npulses=ceil((obj.endTime-obj.startTime)*obj.stimFrequency);
            %Variance of the stimulation timing
            if stimFrequencySD>=0 && stimFrequencySD<100
                obj.stimFrequencySD=stimFrequencySD;
                timeMarkers=(obj.stimFrequencySD*randn(npulses,1)+1/obj.stimFrequency)/obj.h;
                obj.stimTimeMarkers=cumsum(timeMarkers)+obj.startTime/obj.h;
            elseif stimFrequencySD==100 %Poisson process
				obj.stimFrequencySD=stimFrequencySD;
                timeMarkers=(exprnd(1/obj.stimFrequency,npulses,1))/obj.h;
                obj.stimTimeMarkers=cumsum(timeMarkers)+obj.startTime/obj.h;
                if obj.stimTimeMarkers(end)>obj.endTime/obj.h
					obj.stimTimeMarkers(end)=(obj.endTime-0.1)/obj.h;
                end
            elseif stimFrequencySD>1000
                %Phase-locked stimulation with central frequency
                %then stimFrequencySD: number of pulses
                obj.stimTimeMarkers=stimFrequencySD-1000;
                fprintf('Config: Set quantity of stimulus to %f\n',stimFrequencySD-1000);
            elseif obj.stimFrequency==-1
            %Phase-locked stimulation without central frequency
            %then stimFrequencySD: number of pulses
                if stimFrequencySD>0
                    obj.stimTimeMarkers=stimFrequencySD;
                else
                    fprintf('Config: Set quantity of stimulus to 100');
                    obj.stimTimeMarkers=100;
                end
            else
				%Derterministic stimulation
                obj.stimFrequencySD=0;
                timeMarkers=ones(npulses,1)/obj.stimFrequency;
                obj.stimTimeMarkers=cumsum(timeMarkers)+obj.startTime/obj.h;
            end

            if sigmaE>0 && sigmaI>0
                if sigmaI>sigmaE
                    obj.sigmaE=sigmaE;
                    obj.sigmaI=sigmaI;
                elseif sigmaE==sigmaI
                    fprintf('Config: sigmaI muste greater than sigmaE, added 1\n');
                    obj.sigmaE=sigmaE;
                    obj.sigmaI=sigmaE+1;
                %Uncomment to disallow cental inhibition with lateral excitation;
                %else
                 %   fprintf('sigmaI muste greater than sigmaE, values interchanged\n');
                  %  obj.sigmaE=sigmaI;
                  %  obj.sigmaI=sigmaE;
                end
            else
                fprintf('Config: Gaussain deviation must be defined positive, set default values\n');
                obj.sigmaE=1;
                obj.sigmaI=2;
            end
            if min(stimX)>0 && max(stimX)<obj.Nx
                obj.stimX=stimX;
            else
                obj.stimX=1;
            end

            if min(stimY)>0 && max(stimY)<obj.Ny
                obj.stimY=stimY;
            else
                obj.stimY=1;
            end
            if shapeNoiseColor>=0 && shapeNoiseColor<=2
                obj.shapeNoiseColor=shapeNoiseColor;
                %0 None
                %1 White
                %2 Pink
            else
                fprintf('Config: Set shape noise as: None \n');
                obj.shapeNoiseColor=0;
            end
        end
        %Set time parameters of plasticity
        function obj=setPlasticityParams(obj,ratesA,ratesB,tp)
			size_ic=size(ratesA);
			size_id=size(ratesB);
            if size_ic(1)==8 && size_ic(2)==1
                obj.ratesA=ratesA;
            elseif size_ic(1)==1 && size_ic(2)==8
                obj.ratesA=ratesA';
            else
                fprintf('Config: A 8x1 or 1x8 vector expected \n')
                fprintf('Config: Set default values by Asezadeh 2018 \n')
                obj.ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
            end
            if size_ic(1)==8 && size_ic(2)==1
                obj.ratesB=ratesB;
            elseif size_ic(1)==1 && size_ic(2)==8
                obj.ratesB=ratesB';
            else
                fprintf('Config: A 8x1 or 1x8 vector expected \n')
                fprintf('Config: Set default values by Asezadeh 2018 \n')
                obj.ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45];
            end
            obj.initialChangeRates=zeros(8,1);
            eirs_changes=util.eirsConnections(obj.initialChangeRates);
            obj.plasticityMatrix=util.expandConnections(obj.Nx,obj.Ny,eirs_changes); 
			obj.ratesA=ratesA;
			obj.ratesB=ratesB;
			obj.tp=tp;
		end
    end
end
