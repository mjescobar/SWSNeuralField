classdef Simulator
    %"""Solver class, to define parameters of the model and to tag results"""
    %% Properties
    properties
        finalTime
        previousTime
        h
        samplingFrequency
        dx
        second
        delayTCsamples
        model
        modelAAS
        stimulator
        integrator
        monitor
        plasticity
        bufferLength
        samplingLeap
        storageLeap
        storageRows
        plasticityFlag
        steadyV
        steadyQ
        steadyphi
        steadyVm
        steadyVq
        steadyH
        Ach
        indexMax
        indexesN
        indexesE
        indexesI
        indexesR
        indexesS
        indexesDynamic
        indexMaxDynamic
        indexesNoPropagation
        strengthsRatios
        
    end
    %Methods
    methods
        %% Constructor
        function obj=Simulator_old(config,model,integrator,stimulator,monitor,plasticity)
            obj.previousTime=config.previousTime;
            obj.finalTime=config.finalTime;
            obj.h=config.h;
            obj.samplingFrequency=config.samplingFrequency;
            obj.plasticityFlag=config.plasticity_onoff;
            obj.dx=config.dx;            
            obj.samplingLeap=config.samplingLeap;
            obj.storageLeap=config.storageLeap;
            obj.storageRows=config.storageRows;
            obj.model=model;
            obj.stimulator=stimulator; %Stimulation 
            obj.integrator=integrator;
            obj.monitor=monitor; %Store and show data
            obj.plasticity=plasticity; %Dynamical changes of the model
            obj.bufferLength=ceil(1/config.h);
            obj.indexMax=(config.M)*config.Nx*config.Ny;
            obj.indexMaxDynamic=(config.M-1)*config.Nx*config.Ny;
            obj.indexesDynamic=1:(config.M-1)*config.Nx*config.Ny;
            obj.indexesN=(config.M-1)*config.Nx*config.Ny+1:obj.indexMax;
            obj.indexesE=1:config.Nx*config.Ny;
            obj.indexesI=config.Nx*config.Ny+1:2*config.Nx*config.Ny;
            obj.indexesR=2*config.Nx*config.Ny+1:3*config.Nx*config.Ny;
            obj.indexesS=3*config.Nx*config.Ny+1:4*config.Nx*config.Ny;
            obj.indexesNoPropagation=config.Nx*config.Ny+1:(config.M-1)*config.Nx*config.Ny;
            obj.steadyV=0;
            obj.steadyQ=0;
            obj.steadyphi=0;
            obj.strengthsRatios=zeros(8,1);
            %obj.steadyVm=config.steadyVm
            %obj.steadyVq=config.steadyVq
			%obj.steadyH=config.steadyH
            obj.delayTCsamples=ceil(obj.model.t0/(obj.h*2));
            fprintf('Simulator: Instatiation completed \n');
        end
        
        %% Constructor
        function obj=Simulator(config,model,modelAAS,integrator,stimulator,monitor,plasticity)
            obj.previousTime=config.previousTime;
            obj.finalTime=config.finalTime;
            obj.h=config.h;
            obj.samplingFrequency=config.samplingFrequency;
            obj.plasticityFlag=config.plasticity_onoff;
            obj.dx=config.dx;            
            obj.samplingLeap=config.samplingLeap;
            obj.storageLeap=config.storageLeap;
            obj.storageRows=config.storageRows;
            obj.model=model;
            obj.modelAAS=modelAAS;
            obj.stimulator=stimulator; %Stimulation 
            obj.integrator=integrator;
            obj.monitor=monitor; %Store and show data
            obj.plasticity=plasticity; %Dynamical changes of the model
            obj.bufferLength=ceil(1/config.h);
            obj.indexMax=(config.M)*config.Nx*config.Ny;
            obj.indexMaxDynamic=(config.M-1)*config.Nx*config.Ny;
            obj.indexesDynamic=1:(config.M-1)*config.Nx*config.Ny;
            obj.indexesN=(config.M-1)*config.Nx*config.Ny+1:obj.indexMax;
            obj.indexesE=1:config.Nx*config.Ny;
            obj.indexesI=config.Nx*config.Ny+1:2*config.Nx*config.Ny;
            obj.indexesR=2*config.Nx*config.Ny+1:3*config.Nx*config.Ny;
            obj.indexesS=3*config.Nx*config.Ny+1:4*config.Nx*config.Ny;
            obj.indexesNoPropagation=config.Nx*config.Ny+1:(config.M-1)*config.Nx*config.Ny;
            obj.steadyV=0;
            obj.steadyQ=0;
            obj.steadyphi=0;
            obj.Ach=config.Ach;
            obj.strengthsRatios=zeros(8,1);
            %obj.steadyVm=config.steadyVm
            %obj.steadyVq=config.steadyVq
			%obj.steadyH=config.steadyH
            obj.delayTCsamples=ceil(obj.model.t0/(obj.h*2));
            fprintf('Simulator: Instatiation completed \n');
        end
        
        
        %%%% Solve method %%%%%%%%%%%
        %% Principal encapsulating method
        function obj=solve(obj,seed,V,Q,phi,previousV,currentQ)
			%Stimulation control flags
			stimFlags=zeros(1,14);
			noiseFlags=zeros(1,3);
			stimFlags(2)=1;
			%Store flag
            storeFlag=0;
            nstore=1;
            %Storage matrices
            Qstorage=zeros(obj.indexMaxDynamic,obj.storageRows);
            %Potentials vector
            Phistorage=zeros(obj.indexMax,obj.storageRows);
            Vstorage=zeros(obj.indexMaxDynamic,obj.storageRows);
            Phasestorage=zeros(3,obj.storageRows);
            stimstorage=zeros(size(obj.indexesE,2),obj.storageRows);
            markerstorage=zeros(1,obj.storageRows);
            currentTimestorage=zeros(1,obj.storageRows);
            %Random number generator seed
            rng(seed)
            %Mean spatial value of Phi_e, required for closed-loop
            obj.stimulator=obj.stimulator.setMeanPhi(mean(phi(obj.indexesE,1)));
			%All time
			
			for n=1:ceil(obj.finalTime*1/obj.h)
				%% Simulation
				%Set indexes for current step
                [current,currentQ,currentPhi,delay,nstore,samplingFlag,storeFlag]=obj.indexShifting(n,currentQ,nstore);
				%Generate Noise
				[backgroundNoise,noiseFlags]=obj.stimulator.generateNoise(noiseFlags);
				%Generate the stim
				[obj.stimulator,stim,marker,stimFlags]=obj.stimulator.generateStim(n,stimFlags);
				%Set the input
			    phi(obj.indexesN,current)=stim+backgroundNoise;
			    %Perform the simulation step
			    [V,Q,phi,previousV]=obj.step(V,Q,phi,previousV,current,currentQ,currentPhi,delay);
			    %Time in seconds
                currentTime=n*obj.h; %in samples
                %Calculate the spatial mean of phi_e
                stimFlags(9)=mean(phi(obj.indexesE,current));
                
				%Save data
				if samplingFlag==1
					%saveTime=tic();
					Vstorage(:,nstore)=V;
					Qstorage(:,nstore)=Q(:,currentQ(1));
					Phistorage(:,nstore)=phi(:,current);
					Phasestorage(:,nstore)=[obj.stimulator.x_now(1),obj.stimulator.previous_envelope,obj.stimulator.phase];
					stimstorage(:,nstore)=stim;
					markerstorage(nstore)=marker;
					currentTimestorage(nstore)=currentTime;
					%toc(saveTime);
					if storeFlag==1
						obj.monitor.savebatch(Vstorage,Qstorage,Phistorage,stimstorage,markerstorage,currentTimestorage);
						obj.monitor.savePhasebatch(Phasestorage);
						
						if obj.plasticityFlag==2
							[obj.strengthsRatios,meanPhi]=obj.plasticity.update(Phistorage);
							obj.monitor.saveStrengths(obj.model.connectionVector,meanPhi,obj.strengthRatios)
						end
					end
				end
				
				%Plasticity
				if obj.plasticityFlag==2
						obj.model=obj.model.updateConnections(1,obj.storageLeap,obj.strengthsRatios);
				end
				
                 %%Messages
                 if mod(currentTime,5)==0
                    fprintf('Transcurred: %.2f seconds of %s \n',currentTime,obj.monitor.filename);
                    fprintf('Phi_e %.6e  V_e %.6e \n',phi(1,current),V(1,1));
                 end
			end	
		end
		
		%%%% Solve method with AAS%%%%%%%%%%%
        %% Principal encapsulating method
        function obj=solveAAS(obj,seed,V,Q,phi,previousV,Vaas,Qaas,currentQ)
			%Stimulation control flags
			stimFlags=zeros(1,14);
			noiseFlags=zeros(1,3);
			stimFlags(2)=1;
			%Store flag
            storeFlag=0;
            nstore=1;
            %Storage matrices
            Qstorage=zeros(obj.indexMaxDynamic,obj.storageRows);
            Phistorage=zeros(obj.indexMax,obj.storageRows);
            Vstorage=zeros(obj.indexMaxDynamic,obj.storageRows);
            Phasestorage=zeros(3,obj.storageRows);
            stimstorage=zeros(size(obj.indexesE,2),obj.storageRows);
            markerstorage=zeros(1,obj.storageRows);
            currentTimestorage=zeros(1,obj.storageRows);
            %Random number generator seed
            rng(seed)
            %Mean spatial value of Phi_e, required for closed-loop
            obj.stimulator=obj.stimulator.setMeanPhi(mean(phi(obj.indexesE,1)));
			%Plasticiy update(refresh) times
			update_counter=0;
			plasticityUpdateLeap=500;
			%Start simulation
			for n=1:ceil(obj.finalTime*1/obj.h)
				%% Simulation
				%Set indexes for current step
                [current,currentQ,currentPhi,delay,nstore,samplingFlag,storeFlag]=obj.indexShifting(n,currentQ,nstore);
				%Generate Noise
				[backgroundNoise,noiseFlags]=obj.stimulator.generateNoise(noiseFlags);
				%Generate the stim
				[obj.stimulator,stim,marker,stimFlags]=obj.stimulator.generateStim(n,stimFlags);
				%Set the input
			    phi(obj.indexesN,current)=stim+backgroundNoise;
			    %Perform the simulation step
			    %Time in seconds
                currentTime=n*obj.h; %in samples
			    [V,Q,phi,previousV,Vaas,Qaas]=obj.stepAAS(V,Q,phi,previousV,Vaas,Qaas,current,currentQ,currentPhi,delay,currentTime);
			    
                %Calculate the spatial mean of phi_e
                stimFlags(9)=mean(phi(obj.indexesE,current));
                
				%Save data
				if samplingFlag==1
					Vstorage(:,nstore)=V;
					Qstorage(:,nstore)=Q(:,currentQ(1));
					Phistorage(:,nstore)=phi(:,current);
					stimstorage(:,nstore)=stim;
					markerstorage(nstore)=marker;
					currentTimestorage(nstore)=currentTime;
					if storeFlag==1
						obj.monitor.savebatch(Vstorage,Qstorage,Phistorage,stimstorage,markerstorage,currentTimestorage);
						obj.monitor.savePhasebatch(Phasestorage);
						if obj.plasticityFlag==2
							[obj.strengthsRatios,meanPhi]=obj.plasticity.update(Phistorage);
							obj.monitor.saveStrengths(obj.model.connectionVector,meanPhi,obj.strengthsRatios)
						end
					end
				end
				
				%Plasticity
				if obj.plasticityFlag==2
						update_counter=update_counter+1;
						if update_counter>=plasticityUpdateLeap
							obj.model=obj.model.updateConnections(plasticityUpdateLeap,obj.storageLeap,obj.strengthsRatios);
							update_counter=0;
						end
				end
				
				
                 %%Messages
                 if mod(currentTime,5)==0
                    fprintf('Transcurred: %.2f seconds of %s \n',currentTime,obj.monitor.filename);
                    fprintf('Phi_e %.6e  V_e %.6e \n',phi(1,current),V(1,1));
                 end
			end
			
		end
		
		
		
        %% Principal encapsulating method including presolve
        function obj=solveWithPresolve(obj,seed)
			stimFlags=zeros(1,14);
			stimFlags(2)=1;
			noiseFlags=zeros(1,3);
			plasticityNewValuesFlag=0;
			plasticityMarker=0;
			presolverTime=tic();
            storeFlag=0;
            nstore=1;
            %Storage matrices
            Qstorage=zeros(obj.indexMaxDynamic,obj.storageRows);
            Phistorage=zeros(obj.indexMax,obj.storageRows);
            Vstorage=zeros(obj.indexMaxDynamic,obj.storageRows);
            Phasestorage=zeros(3,obj.storageRows);
            stimstorage=zeros(size(obj.indexesE,2),obj.storageRows);
            markerstorage=zeros(1,obj.storageRows);
            currentTimestorage=zeros(1,obj.storageRows);
			[V,Q,phi,previousV,current,currentQ,currentPhi,delay]=obj.presolve(obj.stimulator,seed);
			obj.stimulator=obj.stimulator.setMeanPhi(mean(phi(obj.indexesE,1)));
            fprintf('Presolver completed \n')
			toc(presolverTime);
			for n=1:ceil(obj.finalTime*1/obj.h)
				%Prepare indexes for next step
                [current,currentQ,currentPhi,delay,nstore,samplingFlag,storeFlag]=obj.indexShifting(n,currentQ,nstore);
			    %Input noise
			    [backgroundNoise,noiseFlags]=obj.stimulator.generateNoise(noiseFlags);
				%Stimulus
				[obj.stimulator,stim,marker,stimFlags]=obj.stimulator.generateStim(n,stimFlags);
			    phi(obj.indexesN,current)=stim+backgroundNoise;
			    %Model step
			    [V,Q,phi,previousV]=obj.step(V,Q,phi,previousV,current,currentQ,currentPhi,delay);
                currentTime=n*obj.h; %in samples
                stimFlags(9)=mean(phi(obj.indexesE,current));
				%Save data
			
				if samplingFlag==1
					Vstorage(:,nstore)=V;
					Qstorage(:,nstore)=Q(:,currentQ(1));
					Phistorage(:,nstore)=phi(:,current);
					stimstorage(:,nstore)=stim;
					markerstorage(nstore)=marker;
					currentTimestorage(nstore)=currentTime;
					if storeFlag==1
						obj.monitor.savebatch(Vstorage,Qstorage,Phistorage,stimstorage,markerstorage,currentTimestorage);
						obj.monitor.savePhasebatch(Phasestorage);
						if obj.plasticityFlag==2
							[obj.strengthsRatios,meanPhi]=obj.plasticity.update(Phistorage);
							obj.monitor.saveStrengths(obj.model.connectionVector,meanPhi,obj.strengthsRatios)
						end
					end
				end

				
				%TODO: closed-loop control
				 %stimFlags=obj.stimulator.phaseLockSO(V);
                 
                 %%Messages
                 if mod(currentTime,5)==0
                    fprintf('Transcurred: %.2f seconds of %s \n',currentTime,obj.monitor.filename);
                    fprintf('Phi_e %.6e  V_e %.6e \n',phi(1,current),V(1,1));
                 end 
			end	%Simulation step
		end %function 
		
	%%%%%%%%%%%%%%%%%%%%
	%% Presolver method, 
        % Start near the steady state points in all
		%populations
		function[steadyV,steadyQ,steadyPhi,steadyPreviousV,current,currentQ,currentPhi,delay]=presolve(obj,stimulator,seed)
			timeSkip=ceil(0.5*1/obj.h);
            timeMean=ceil((obj.previousTime-0.5)*1/obj.h);
            %%% Execute the solver for time_presolver seconds before save results
            %Firing response vector
            Q=zeros(obj.indexMaxDynamic,3);
            %Potentials vector
            phi=zeros(obj.indexMax,obj.bufferLength);
            V=zeros(obj.indexMaxDynamic,1);
            previousV=zeros(obj.indexMaxDynamic,2);
            currentQ=[3,1,2];
            rng(seed);            
            backgroundNoise=stimulator.noiseMean+stimulator.noiseSD*randn(obj.indexMax/5,ceil(1/obj.h)*obj.previousTime);
            meanV=0;
            meanphi=0;
            minv=1;
            meanQ=0;
            SOth=-4e-5;
            nstore=0;
            %Simulation for the previousTime
            for n=1:ceil(1/obj.h*obj.previousTime)
				[current,currentQ,currentPhi,delay,nstore,samplingFlag,storeFlag]=obj.indexShifting(n,currentQ,nstore);
				phi(obj.indexesN,current)=backgroundNoise(:,n);
				%Simulation step
                [V,Q,phi,previousV]=obj.step(V,Q,phi,previousV,current,currentQ,currentPhi,delay);
                %Save values and search for the minimum value
                if n>timeSkip
					%Means of coritcal population
					meanV=meanV+mean(V(obj.indexesE,1));
					meanQ=meanQ+mean(Q(obj.indexesE,currentQ(1)));
					meanphi=meanphi+mean(phi(obj.indexesE,current));
					minv_candidate=min(V(obj.indexesE,1));
					if minv_candidate<minv
						minv=minv_candidate;
					end
					if minv<-1e-4
					   SOth=-4e-5;
					else
					   SOth=0.75*(minv-meanV);
					end
                end
            end %end for
            meanV=meanV/timeMean;
            meanQ=meanQ/timeMean;
            meanphi=meanphi/timeMean;
            obj=setSteadyValues(obj,meanV,meanQ,meanphi);
            obj.stimulator.setSOth(SOth);
            steadyV=V;
            steadyPhi=phi;
            steadyQ=Q;
            steadyPreviousV=previousV;
		end
		%% Presolver method with AAS system
        % Start near the steady state points in all
		%populations
		function[steadyV,steadyQ,steadyPhi,steadyPreviousV,Vaas,Qaas,current,currentQ,currentPhi,delay]=presolveAAS(obj,stimulator,seed)
			timeSkip=ceil(0.5*1/obj.h);
            timeMean=ceil((obj.previousTime-0.5)*1/obj.h);
            %%% Execute the solver for time_presolver seconds before save results
            %Firing response vector
            Q=zeros(obj.indexMaxDynamic,3);
            Vaas=[-0.54e-3,-3.75e-3,13.65e-9];
            Qaas=obj.modelAAS.firingResponseAAS(Vaas);
            %Potentials vector
            phi=zeros(obj.indexMax,obj.bufferLength);
            V=zeros(obj.indexMaxDynamic,1);
            previousV=zeros(obj.indexMaxDynamic,2);
            currentQ=[3,1,2];
            rng(seed);            
            backgroundNoise=stimulator.noiseMean+stimulator.noiseSD*randn(obj.indexMax/5,ceil(1/obj.h)*obj.previousTime);
            meanV=0;
            meanphi=0;
            minv=1;
            meanQ=0;
            SOth=-4e-5;
            nstore=1;
            %flags
            nstore=0;
            %Simulation for the previousTime
            for n=1:ceil(1/obj.h*obj.previousTime)
				[current,currentQ,currentPhi,delay,nstore,samplingFlag,storeFlag]=obj.indexShifting(n,currentQ,nstore);
				phi(obj.indexesN,current)=backgroundNoise(:,n);
				%Simulation step
                [V,Q,phi,previousV,Vaas,Qaas]=obj.stepAAS(V,Q,phi,previousV,Vaas,Qaas,current,currentQ,currentPhi,delay,n*obj.h);
                %Save values and search for the minimum value
                if n>timeSkip
					%Means of coritcal population
					meanV=meanV+mean(V(obj.indexesE,1));
					meanQ=meanQ+mean(Q(obj.indexesE,currentQ(1)));
					meanphi=meanphi+mean(phi(obj.indexesE,current));
					minv_candidate=min(V(obj.indexesE,1));
					if minv_candidate<minv
						minv=minv_candidate;
					end
					if minv<-1e-4
					   SOth=-4e-5;
					else
					   SOth=0.75*(minv-meanV);
					end
                end
            end %end for
            meanV=meanV/timeMean;
            meanQ=meanQ/timeMean;
            meanphi=meanphi/timeMean;
            obj=setSteadyValues(obj,meanV,meanQ,meanphi);
            obj.stimulator.setSOth(SOth);
            steadyV=V;
            steadyPhi=phi;
            steadyQ=Q;
            steadyPreviousV=previousV;
		end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	%% Step
        %kernel of the simulations
		function [V,Q,phi,previousV]=step(obj,V,Q,phi,previousV,current,currentQ,currentPhi,delay)
			%Dendrytic Voltages, taking in account the thalamocortical-delay
			%fprintf('t1:\n');
			%t1=tic();             
			%vAux=obj.model.connectionMatrixNodelay*phi(:,current)+obj.model.connectionMatrixDelay*phi(:,delay);
			vAux=prodsum(phi(:,current),obj.model.connectionMatrixNodelay,phi(:,delay),obj.model.connectionMatrixDelay);
            %toc(t1);
			%Temporal dynamics of dendrytic Voltages
			%t2=tic();
			previousV=obj.integrator.stepEuler(previousV,vAux(obj.indexesDynamic));
			%fprintf('t2:\n');
			%toc(t2);
			%store the n+1 values for the n+1 step
			%Store soma voltages (squeeze(vab))
            V(:,1)=previousV(:,1);
            %t3=tic();
            Q(:,currentQ(1))=obj.model.firingResponse(V);
            %fprintf('t3:\n');
            %toc(t3);
            %Spatial propagation Q:1:3 Phi 2:3
            %t4=tic();
            phi(obj.indexesE,currentPhi(1))=obj.integrator.propagator.propagate(Q(obj.indexesE,currentQ),phi(obj.indexesE,currentPhi(2:3)));
            phi(obj.indexesNoPropagation,currentPhi(1))=Q(obj.indexesNoPropagation,currentQ(1));
            %fprintf('Phie:, %.8f,  Phin:, %.8f\n',phi(1,current),phi(1025,current)); 
            %fprintf('n:%d n+1:%d \n', current, currentPhi(1));
            %fprintf('post Phi(%d):%.9e; Phi(%d):%.9e\n', currentPhi(1), phi(1,currentPhi(1)), currentPhi(2), phi(1,currentPhi(2)));
            %fprintf('t4:\n');
            %toc(t4);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %StepAAS
        %kernel of the simulations
		function [V,Q,phi,previousV,Vaas,Qaas]=stepAAS(obj,V,Q,phi,previousV,Vaas,Qaas,current,currentQ,currentPhi,delay,time)
			%Dendrytic Voltages, taking in account the thalamocortical-delay
			phi_current=phi(:,current)+Qaas(2)*obj.model.stencilAAS;       
			%vAux=obj.model.connectionMatrixNodelay*phi_current+obj.model.connectionMatrixDelay*phi(:,delay);
			vAux=prodsum(phi_current,obj.model.connectionMatrixNodelay,phi(:,delay),obj.model.connectionMatrixDelay);
			D=obj.modelAAS.sleepDrive(time,Vaas(3));
			u=obj.modelAAS.vector_u(Qaas,D,obj.Ach);
			
			%Soma voltages
			[previousV,Vaas]=obj.integrator.stepEulerAAS(previousV,vAux(obj.indexesDynamic),Vaas,u);
			
			%Store soma voltages 
            V(:,1)=previousV(:,1);
            %Firing response
            Q(:,currentQ(1))=obj.model.firingResponse(V);
            Qaas=obj.modelAAS.firingResponseAAS(Vaas);
            %Spatial Propagation
            phi(obj.indexesE,currentPhi(1))=obj.integrator.propagator.propagate(Q(obj.indexesE,currentQ),phi(obj.indexesE,currentPhi(2:3)));
            phi(obj.indexesNoPropagation,currentPhi(1))=Q(obj.indexesNoPropagation,currentQ(1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=setSteadyValues(obj,meanV,meanQ,meanphi)
			%% Set values
            obj.steadyV=meanV;
            obj.steadyQ=meanQ;
            obj.steadyphi=meanphi;
        end
  
  
  	%%%%%%%%%%%%%%%%%%%%%%%%%%
  	%Control of timing indexes
		function [current,currentQ,currentPhi,delay,nstore,samplingFlag,storageFlag]=indexShifting(obj,n,currentQ,nstore)
			%%Return the current position to store temporally the data, the delay sample
			%and a flag to store permanently the data
			%
			current=mod((n-1+obj.delayTCsamples),obj.bufferLength)+1;
			if current>obj.delayTCsamples
				delay=current-obj.delayTCsamples;
			else
				delay=obj.bufferLength-obj.delayTCsamples+current;
			end
			%Indexes for propagation (after the integration step but in same solver step)[n+1, n, n-1];
            
            if current==1
                currentPhi=[2, 1, obj.bufferLength];  
            elseif current<obj.bufferLength && current>1
                currentPhi=[current+1, current, current-1];
            else
                currentPhi=[1, obj.bufferLength, obj.bufferLength-1];
            end
            %indexes of Q: 123->231->312->123
            %first index the n+1 samples Q(1)=Qn+1, Q(2)=Qn+2, Q(3)=Qn+3,
            %at the third step Q([3,1,2])=[Qn+3, Qn+2, Qn+1], the desired order.
            currentQ=circshift(currentQ,1);
			%Monitor or store data timing control
			samplingFlag=0;
			storageFlag=0;
			sample=mod(n,obj.samplingLeap);
			store=mod(n,obj.storageLeap);
			if sample==0
				samplingFlag=1;
				nstore=nstore+1;
				if store==0
					%Plasticity and storage timing control
					storageFlag=1;
					nstore=1;
				end
			end
			
			
		end
    end
    
end

