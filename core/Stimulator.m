classdef Stimulator
    
    %Properties of the stimulator
    properties
        stimMode
        stimFrequency %
        stimFrequencySD %10%
        stimTimeMarkers
        stimPulseDuration
        stimPulseDurationSamples
        samplingLeap
        startTime
        endTime
        dt
        stimAmplitude
        stimAmplitudeSD
        noiseMean
        noiseSD
        noiseColor
        shapeNoiseColor
        targetPhase
        stimShape %First test rectangular
        Nx
        Ny
        sigmaE
        sigmaI
        sigmaG
        stimX
        stimY
        SOth
        numDelta
        denDelta
        numEnvelope
        denEnvelope
        numShift
        denShift
        previous_envelope_filter
        previous_delta_filter
        previous_shift_filter
        previous_alpha_power_filter
        previous_theta_power_filter
        previous_beta_power_filter
        previous_eeg_power_filter
        previous_delta
        previous_envelope
        previous_alpha_power
        previous_theta_power
        previous_beta_power
        previous_eeg_power
        x_now
        y_now
        origen_eeg
        shiftphase_eeg
        phase
        ENVELOPE_OFFSET
        COUNT_STIMULI
        gain_shiftphase
        phaseOffset
        timePauseStart
        pauseTime
        meanPhi
        npulses
        interpulsetime
        peak_value
        envelope_accum
        envelope_count
        envelope_flag
        peak_value_shift
        envelope_accum_shift
        envelope_count_shift
        envelope_flag_shift
        previous_envelope_shift
        envelope_tol
        envelope_tol_shift
        envelope_wait
        envelope_wait_shift
    end
    methods
        %Constructor
        function obj=Stimulator(mode,config)
            %%%These come from previous defined parameters
            %dt: h, solver step
            %Nx,Ny,
            %%%These should be constants defined at the start
            %sigmaE: excitatory standard deviation
            %sigmaI: inhibitory standard deviation
            %sigma1>sigma2: central inhibitory,
            %lateral excitatory, sigma2>sigma1: central excitatory, lateral
            %inhibitory.
            %theta: stimuli threshold
            %x,y: position in the grid Nx*Ny
            %duration_stim time ON of stimuli (samples)
            %frequency: frequency of stimulation
            %amplitude: stimulus amplitude
            
            %Mode 0: open loop, 1: closed loop, 2: closed loop-modification, 3:sham
            if mode<4
                obj.stimMode=mode;
            else
                obj.stimMode=0;
            end
            
            obj.npulses=config.npulses;
            obj.interpulsetime=config.interpulsetime;
            obj.Nx=config.Nx;
            obj.Ny=config.Ny;
            obj.dt=config.h;
            obj.stimFrequency=config.stimFrequency;
            obj.stimFrequencySD=config.stimFrequencySD;
            obj.stimTimeMarkers=config.stimTimeMarkers;
            obj.stimPulseDuration=config.stimPulseDuration;
            obj.stimPulseDurationSamples=config.stimPulseDurationSamples;
            obj.samplingLeap=config.samplingLeap;
            obj.startTime=config.startTime;
            obj.endTime=config.endTime;
            obj.stimAmplitude=config.stimAmplitude;
            obj.stimAmplitudeSD=config.stimAmplitudeSD;
            obj.noiseMean=config.noiseMean;
            obj.noiseSD=config.noiseSD;
            obj.noiseColor=config.noiseColor; %1 White noise, 2 Pink noise
            obj.shapeNoiseColor=config.shapeNoiseColor;
            obj.targetPhase=config.targetPhase;
            obj.sigmaE=config.sigmaE;
            obj.sigmaI=config.sigmaI;
            obj.stimX=config.stimX;
            obj.stimY=config.stimY;
            obj.phase=0;
            %obj.SOth=-4e-5;
            obj.SOth=-0.1;
            obj.sigmaG=(obj.stimPulseDurationSamples/2)^2;
            %shape: 0: rectangular, 1: rise ramp, 2: decrease ramp, 3: triangle,
            %4: gaussian, 5: Hamming, 6:sinc, 7:sin, 8:sin2,  9:cos, 10: exponential long-tail,
            %11: exponential short_tail
            obj.stimShape=config.stimShape;
            
            [obj.numDelta, obj.denDelta,obj.numShift,obj.denShift]=util.filterCoeffs(obj.stimFrequency,1/config.h);
            %if config.h==1e-4
			%	fprintf('Stimulator: Filters coeficients for h=1e-04\n')
			%	obj.numDelta=[1e-4, 1e-4, -1e-4, -1e-4];
			%	obj.denDelta=[1.0, -1.999799605253469, 0.99980001, 0];
			%	%obj.numEnvelope=[2.46685308301412e-08, 4.93370616602824e-08, 2.46685308301412e-08];
			%	%obj.denEnvelope=[1, -1.99955571171349, 0.999555810387613];
			%	obj.numShift=[1.0,  1.056428578934086,   0.390785721040677,   0.225]*7.07e-11;
			%	obj.denShift=[1.0, -2.998799605253469,  2.997599815648216, -0.99880020999];
			%end
			delta_order=size(obj.numDelta);
			delta_order=delta_order(2);
			obj.previous_delta_filter=zeros(delta_order-1,1);
			obj.previous_delta=zeros(delta_order,1);
			obj.previous_shift_filter=zeros(delta_order-1,1);
			obj.shiftphase_eeg=zeros(delta_order,1);
			obj.previous_envelope=0;
			obj.envelope_count=0;
			obj.envelope_accum=0;
			obj.envelope_flag=0;
			obj.envelope_tol=5e-4;
            obj.envelope_wait=false;
			obj.previous_envelope_shift=0;
			obj.envelope_count_shift=0;
			obj.envelope_accum_shift=0;
			obj.envelope_flag_shift=0;
			obj.envelope_tol_shift=5e-4;
            obj.envelope_wait_shift=false;
			obj.previous_theta_power_filter=zeros(delta_order-1,1);
			obj.previous_theta_power=zeros(delta_order,1);
			obj.previous_alpha_power_filter=zeros(delta_order-1,1);
			obj.previous_alpha_power=zeros(delta_order,1);
			obj.previous_beta_power_filter=zeros(delta_order-1,1);
			obj.previous_beta_power=zeros(delta_order,1);
			obj.previous_eeg_power_filter=zeros(delta_order-1,1);
			obj.previous_eeg_power=zeros(delta_order,1);
			obj.x_now=zeros(3,1);
			obj.y_now=zeros(3,1);
			obj.ENVELOPE_OFFSET=4e-3;
			obj.peak_value=obj.ENVELOPE_OFFSET;
			obj.peak_value_shift=obj.ENVELOPE_OFFSET;
			obj.COUNT_STIMULI=0;
			obj.origen_eeg=zeros(3,1);
			obj.shiftphase_eeg=zeros(delta_order,1);
			obj.gain_shiftphase=config.samplingFrequency/(2*pi);
			obj.phaseOffset=0;
			obj.pauseTime=2; %2 seconds of pause between pulses
			obj.timePauseStart=0;
			obj.meanPhi=12; %approx
			fprintf('Stimulator: Instatiation completed\n')
        end
        
        function obj=setMeanPhi(obj,newMeanPhi)
			obj.meanPhi=newMeanPhi;
		end
        
        function [noise,flags]=generateNoise(obj,flags)	
            noise=zeros(obj.Nx*obj.Ny,1);
			if obj.noiseColor==0
				flags=flags;
            elseif obj.noiseColor==1
                noise=obj.noiseMean+obj.noiseSD*randn(obj.Nx*obj.Ny,1);
            elseif  obj.noiseColor==2
                [noise_node,randgen]=util.pinknoise(flags);
                flags=randgen;
                for k=1:obj.Nx*obj.Ny
                    noise_out=obj.noiseMean+obj.noiseSD*noise_node*rand(1,1);
                    noise(k)=noise_out;
                end
            elseif obj.noiseColor==3
				[noise_node,randgen]=util.pinknoise(flags);
                flags=randgen;
                noise=obj.noiseMean+obj.noiseSD*noise_node*randn(256,1);
            end
        end
        
        function [noise,flags]=generateHomeostaticNoise(obj,n,flags,homeostaticPhase)
            noise=zeros(obj.Nx*obj.Ny,1);
            hMean=pi/3+pi/3*cos(2*pi*(1/5400)*obj.dt*n+homeostaticPhase);
            if obj.noiseColor==1
                noise=hMean+obj.noiseSD*randn(obj.Nx*obj.Ny,1);
            elseif  obj.noiseColor==2
                [noise_node,randgen]=util.pinknoise(flags);
                flags=randgen;
                for k=1:obj.Nx*obj.Ny
                    noise_out=hMean+obj.noiseSD*noise_node*rand(1,1);
                    noise(k)=noise_out;
                end
            end 
        end
        
        function [input]=generateRampInput(obj,n,slope)
            mone=ones(obj.Nx*obj.Ny,1);
            input=n*slope*mone;    
        end
        
        function [shape_now,now_flags]=generateShape(obj,now_flags,n)
			if obj.stimShape==0 %rectagular
				shape_now=obj.stimAmplitude;
			elseif obj.stimShape==1 %rising ramp
				shape_now=obj.stimAmplitude*(n-now_flags(1))/(obj.stimPulseDurationSamples);
			elseif obj.stimShape==2 %decreasing ramp
				shape_now=obj.stimAmplitude-obj.stimAmplitude/(obj.stimPulseDurationSamples)*(n-now_flags(1));
			elseif obj.stimShape==3 %triangular
					shape_now=obj.stimAmplitude-2/obj.stimPulseDurationSamples*obj.stimAmplitude*abs(n-now_flags(1)-obj.stimPulseDurationSamples/2);
			elseif obj.stimShape==4 %Gaussian
				shape_now=obj.stimAmplitude*exp(-(n-now_flags(1)-obj.stimPulseDurationSamples/2)^2/obj.sigmaG);
			elseif obj.stimShape==5 %Hamming
				shape_now=obj.stimAmplitude*(0.53836-0.46164*cos(2*pi*(n-now_flags(1))/(obj.stimPulseDurationSamples-1)));
			elseif obj.stimShape==6 %sinc
				if (n-now_flags(1)-obj.stimPulseDurationSamples/2)==0
					shape_now=obj.stimAmplitude;
				else
					shape_now=obj.stimAmplitude/2*sin(2*pi*obj.dt/obj.stimPulseDuration*(n-now_flags(1)-obj.stimPulseDurationSamples/2))/(pi*obj.dt/obj.stimPulseDuration*(n-now_flags(1)-obj.stimPulseDurationSamples/2));
				end
			elseif obj.stimShape==7 %sin
				shape_now=obj.stimAmplitude*abs(sin(pi*obj.dt/obj.stimPulseDuration*(n-now_flags(1))));
			elseif obj.stimShape==8 %sin2
				shape_now=obj.stimAmplitude*abs(sin(pi*obj.dt/obj.stimPulseDuration*(n-now_flags(1))))^2;
			elseif obj.stimShape==9 %cos
				shape_now=obj.stimAmplitude*abs(cos(pi*obj.dt/(2*obj.stimPulseDuration)*(n-now_flags(1))+3/2*pi));
			elseif obj.stimShape==10 %exp long tail
				shape_now=obj.stimAmplitude*exp(-abs(n-now_flags(1)-obj.stimPulseDurationSamples/2)/(obj.stimPulseDurationSamples/9));
			elseif obj.stimShape==11 %exp short tail
				shape_now=obj.stimAmplitude*exp(-abs(n-now_flags(1)-obj.stimPulseDurationSamples/2)/(obj.stimPulseDurationSamples/4));
			elseif obj.stimShape==12 %Trapeze
				if (n-now_flags(1)-3*obj.stimPulseDurationSamples/4)>0
					shape_now=obj.stimAmplitude-4*obj.stimAmplitude/(obj.stimPulseDurationSamples)*(n-now_flags(1)-3*obj.stimPulseDurationSamples/4);
				elseif (n-now_flags(1)-obj.stimPulseDurationSamples/4)<0
					shape_now=obj.stimAmplitude*4*(n-now_flags(1))/(obj.stimPulseDurationSamples);
				else
					shape_now=obj.stimAmplitude;
				end
			elseif obj.stimShape==13 %Rectangular trapezoid
				if (n-now_flags(1))>obj.stimPulseDurationSamples/2
					shape_now=obj.stimAmplitude-2*obj.stimAmplitude/(obj.stimPulseDurationSamples)*(n-now_flags(1)-obj.stimPulseDurationSamples/2);
				else
					shape_now=obj.stimAmplitude;
				end
			elseif obj.stimShape==14 %Exponential
				shape_now=obj.stimAmplitude*exp(-(n-now_flags(1))/(obj.stimPulseDurationSamples/4));
			end
			
			if obj.shapeNoiseColor==1 %White noise
				shape_now=shape_now*(1+obj.stimAmplitudeSD*randn(1,1));
			elseif obj.shapeNoiseColor==2
				[noise,randgen]=util.pinknoise(now_flags(4:6));
				now_flags(4:6)=randgen;
				shape_now=shape_now*(1+obj.stimAmplitudeSD*noise);
			end
        end
        
        function [obj,stimulus,marker,flags]=generateStim(obj,n,flags)
            %Inputs:
            %These are dynamic
            %n: actual sample
            %flags: quantities relative to the stimulation type, maybe someone is actually a flag.
            %Open-loop: flags=[0,0,0,0,0,0]: start n, prev_value stimulus, sample, rangen0,randgen1,randgen2, old_marker
            %Closed-loop: flags=[0,0,0,0,0,0,0]: start n, prev_value stimulus, sample, rangen0,randgen1,randgen2, old_marker, V, flag_shot, flag_pause, flag_hold, flag_search, pulses,
            %Outputs: 
            
            stimulus=zeros(obj.Nx*obj.Ny,1);
            now_flags=flags;
            marker=[];
            old_marker=now_flags(7);
            %Stimulation is overall nodes
            template=zeros(obj.Nx,obj.Ny);
            %Stimulation frequency (interpulse)
            %Control to start stimulation
            if (n*obj.dt)>=obj.startTime && (n*obj.dt)<=obj.endTime
                if obj.stimMode==0
					%Open-loop mode
					if now_flags(2)<=length(obj.stimTimeMarkers)
						if n>=obj.stimTimeMarkers(now_flags(2))
							now_flags(1)=n;
							now_flags(3)=1;
							marker=1;
							now_flags(2)=now_flags(2)+1;
						end
					end
				elseif obj.stimMode==1	
						%closed loop: Num stims mode
                        if now_flags(2)<=obj.stimTimeMarkers
							[obj,shape_now,marker,stimulationFlags]=phaseLock(obj,now_flags,n);
							now_flags=stimulationFlags;
                        else
                            shape_now=0;
                        end 
				elseif obj.stimMode==2	
						%Closed loop: Ngo Mode
						if now_flags(2)<=obj.stimTimeMarkers
							[obj,shape_now,marker,stimulationFlags]=phaseLockSO(obj,now_flags,n);
							now_flags=stimulationFlags;
						else
							shape_now=0;
						end
				else
					shape_now=0; %Incorrect stimulation mode
				end
			else
				shape_now=0; %Out of stimulation time
			end
            %In general, stim(t)=shape*amplitude(t)*ocurrence(t);
			if now_flags(3)<=obj.npulses && now_flags(3)>0
				if (n-now_flags(1))<=obj.stimPulseDurationSamples && (n-now_flags(1))>=0
					[shape_now,now_flags]=generateShape(obj,now_flags,n);				
				elseif (n-now_flags(1))>obj.stimPulseDurationSamples
					fprintf('New small pulse \n')
					shape_now=0;
					%Advance start-time marker for the next pulse
					now_flags(1)=n+ceil(obj.interpulsetime/obj.dt);
					%Next pulse
					now_flags(3)=now_flags(3)+1;
				end
				%Relu with saturation
				%if shape_now<0
				%	shape_now=0;
				%elseif shape_now>1
				%	shape_now=theta;
				%end
			else
				shape_now=0; %Out of pulse duration time
			end
			%Stim MODE 0
			if (n-now_flags(1))>=(obj.stimPulseDurationSamples*obj.npulses+ceil(obj.interpulsetime/obj.dt)*(obj.npulses-1))
				now_flags(3)=0; %End of stimulus
			end
            %Spatial shape, always a mexican hat or a sum of Gaussians
            for i=1:obj.Nx
                for j=1:obj.Ny
                    if length(obj.stimX)==1
                        template(i,j)=template(i,j)+shape_now/(sqrt(2*pi)*obj.sigmaE)*exp(-((i-obj.stimX)^2+(j-obj.stimY)^2)/obj.sigmaE)-shape_now/(sqrt(2*pi)*obj.sigmaI)*exp(-((i-obj.stimX)^2+(j-obj.stimY)^2)/obj.sigmaI);
                    else
                        for pos=1:length(obj.stimX)
                            template(i,j)=template(i,j)+shape_now/(sqrt(2*pi)*obj.sigmaE)*exp(-((i-obj.stimX(pos))^2+(j-obj.stimY(pos))^2)/obj.sigmaE)-shape_now/(sqrt(2*pi)*obj.sigmaI)*exp(-((i-obj.stimX(pos))^2+(j-obj.stimY(pos))^2)/obj.sigmaI);
                        end
                    end
                end
            end
            for k=1:obj.Nx*obj.Ny
                nx=mod(k,obj.Nx);
                if nx==0
                    nx=obj.Nx;
                end
                ny=floor(k/obj.Ny)+1;
                if ny>obj.Ny
                    ny=obj.Ny;
                end
                stimulus(k)=template(nx,ny);
            end
            if isempty(marker)
				marker=0;
			end
			now_flags(7)=marker;
            flags=now_flags;
        end
        
        
        function [obj,shape_now,marker,stimulationFlags]=phaseLockSONgo(obj,flags,n)
			currentSample=flags(1);
			old_marker=flags(7);
			V=flags(9);
			flag_shot=flags(10);
			flag_pause=flags(11);
			flag_hold=flags(12);
			flag_search=flags(13);
			pulses=flags(14);
			EPSILON_ANGLE=0.05*pi;
			marker=[];
			stimulus=zeros(obj.Nx*obj.Ny,1);
			STIM_SHAM=1;
            obj.previous_delta(2:end)=obj.previous_delta(1:end-1);
            obj.previous_delta(1)=V-obj.meanPhi;
            shape_now=0;
            nsamples=2000; %0.4 second para h=1e-4
            %filtering
            [filtro_delta,obj.previous_delta_filter]=util.filtrar(obj.numDelta,obj.denDelta,obj.previous_delta,obj.previous_delta_filter); %filtrar
            obj.x_now(2:3)=obj.x_now(1:2);
            obj.y_now(2:3)=obj.y_now(1:2);
            %Rectifiying
            if filtro_delta>=0
                obj.x_now(1)=filtro_delta;
            else
                obj.x_now(1)=-filtro_delta;
            end
            
			%New envelope detection method
            if obj.envelope_flag==0 && obj.envelope_count>0 && n>nsamples
                %Calculate a new tolerance afte nsamples of simulated data
                obj.envelope_tol=obj.envelope_accum/obj.envelope_count;
                obj.envelope_flag=1;
                fprintf('Tolerance envelope %.4e\n',obj.envelope_tol);
            end
            if obj.x_now(1)<obj.envelope_tol
                %Set the search for peak value to 0.01 if the signal is below the tolerance
                obj.peak_value=obj.ENVELOPE_OFFSET;
            elseif abs(obj.x_now(1)-obj.x_now(2))<=obj.envelope_tol && obj.x_now(1)>=obj.x_now(2)
                %when the signal arrives to a peak
                if obj.x_now(3)<obj.peak_value && obj.x_now(1)>obj.peak_value 
                    %If the point is higher than the previous value, search for the peak value 
                    obj.envelope_wait=true;
                    if n<nsamples
                        %Add the derivative to calculate the average tolerance value for zero derivative
                        obj.envelope_accum=obj.envelope_accum+abs(obj.x_now(1)-obj.x_now(2));
                        obj.envelope_count=obj.envelope_count+1;
                    end
                else
                    %If the point is not the peak, save the value to search for higher values than it
                    obj.peak_value=obj.x_now(1);
                end
            elseif obj.envelope_wait==true && obj.x_now(1)<obj.x_now(2)
                %If the value is lesser than the previous one, the peak was reached, then save that value as the envelope 
                obj.previous_envelope=obj.x_now(2);
                obj.envelope_wait=false;
            end

            envolvente=obj.previous_envelope;
            if envolvente<obj.ENVELOPE_OFFSET %envelope lesser than 0.01
                envolvente=obj.ENVELOPE_OFFSET; %envelope set to 0.01
                %avoid the stimuli delivery until next oscillation
                flag_hold=1;
            end
            %Save the actual envelope
            obj.previous_envelope=envolvente;
            %delta(i)=filtro_delta;
            %Normalización
            eeg_normalizado=filtro_delta/envolvente;
            
            
            
            %Shifting phase -270 degrees (or forwarding 90 degrees)
            obj.shiftphase_eeg(2:end)=obj.shiftphase_eeg(1:end-1);
            obj.shiftphase_eeg(1)=eeg_normalizado;
            [shiftphase,obj.previous_shift_filter]=util.filtrar(obj.numShift,obj.denShift,obj.shiftphase_eeg,obj.previous_shift_filter);
            
            if shiftphase>=0
                obj.y_now(1)=shiftphase;
            else
                obj.y_now(1)=-shiftphase;
            end
            

            %New envelope detection method
            if obj.envelope_flag_shift==0 && obj.envelope_count_shift>0 && n>nsamples
                obj.envelope_tol_shift=obj.envelope_accum_shift/obj.envelope_count_shift;
                obj.envelope_flag_shift=1;
                fprintf('Tolerance envelope shift %.4e\n',obj.envelope_tol_shift);
            end
            if obj.y_now(1)<obj.envelope_tol_shift
                obj.peak_value_shift=obj.ENVELOPE_OFFSET;
            elseif abs(obj.y_now(1)-obj.y_now(2))<=obj.envelope_tol_shift && (obj.y_now(1)-obj.y_now(2))>=0
                if obj.y_now(3)<obj.peak_value_shift && obj.y_now(1)>obj.peak_value_shift
                    obj.envelope_wait_shift=true;   
                    if n<nsamples
                        obj.envelope_accum_shift=obj.envelope_accum_shift+abs(obj.y_now(1)-obj.y_now(2));
                        obj.envelope_count_shift=obj.envelope_count_shift+1;
                    end
                else
                    obj.peak_value_shift=obj.y_now(1);
                end
            elseif obj.envelope_wait_shift==true && obj.y_now(1)<obj.y_now(2)
                obj.previous_envelope_shift=obj.y_now(2);
                obj.envelope_wait_shift=false;
                
            end

            envolvente_shift=obj.previous_envelope_shift;
            if envolvente_shift<0.25%The envelope of the shifted phase signal at leat must be 0.1
                envolvente_shift=0.25; %Set the envelope to 0.1
            end
            %Save the actual envelope
            obj.previous_envelope_shift=envolvente_shift;

            eeg_normalizado_90=shiftphase/envolvente_shift;
              
			%Calculating the current phase with the arctang
			%the shifted-phase signal is a sine
			%the original signal is takin in account as the cosine
			fase=atan(eeg_normalizado_90/(eeg_normalizado+1e-12));
			
			faseraw=fase; %save for comparisons
			
			%Translate from [-pi/2, pi/2] domain to [0, 2pi] range
			if obj.shiftphase_eeg(1)>0 && obj.shiftphase_eeg(2)<0 %cosine positively crossing zero 
				%aumento_fase=pi/2; python
				obj.phaseOffset=pi/2;
			elseif obj.shiftphase_eeg(2)>0 && obj.shiftphase_eeg(1)<0 %cosine neagtively crossing zero
				%aumento_fase=3/2*pi; python
				obj.phaseOffset=pi*3/2;
				% 								elseif shiftphase_eeg(1)>0 && shiftphase_eeg(2)<0 %coseno subida
				% 									%aumento_fase=5/2*pi; python
				% 									aumento_fase=pi/2;
				% 								elseif shiftphase_eeg(2)>0 && shiftphase_eeg(1)<0 %coseno bajada
				% 									%aumento_fase=pi/2; python
				% 									aumento_fase=pi*3/2;
			end
			%sumar el desfase correspondiente segun el
			%cruce por cerostart
			fase=fase+obj.phaseOffset;
			delta=filtro_delta;
			%Si la fase resulta mayor a 2pi mantener en esa cota
			if fase>2*pi
				obj.phase=2*pi;
			elseif fase<0
				obj.phase=0;
			else
				obj.phase=fase;
			end
			
			%fprintf('%.6e\t %.6e\t %.6e\t %.6e \n',filtro_delta,envolvente,eeg_normalizado,obj.phase);  
			if (obj.previous_delta(1)<obj.SOth  && flag_hold==0)
				flag_search=1;
				%fprintf('flag_search set \n');
			end
			
			if (flag_search==1 && obj.phase>=obj.targetPhase-EPSILON_ANGLE && obj.phase<=obj.targetPhase+EPSILON_ANGLE)%%&& filtro_delta>0
				
				if flag_shot==0
				fprintf('flag_shot set \n');
					flag_shot=1;     %Perform a pulse of the stimulus
					flag_search=0;  %Deactivate for not search the angle again until the first pulse is delivered
				end
			end
			if flag_shot==1
				flag_shot=0;
				flag_hold=1; %Stim delivered this time, wait until next cyvle
				if STIM_SHAM==1
					%If STIM radio button is selected
					if flag_pause==1
						fprintf('On pause \n');
						%Similar to SHAM condition
						if pulses==0
							pulses=1;
							marker=2;
						else
							pulses=3;
							marker=4;
						end
						
					else
						%Stimulation
						fprintf('Stimulus at: %f \n',obj.stimTimeMarkers(flags(2)));
						fprintf('Stim number %d \n',flags(2));
						flags(1)=n;
						flags(3)=1;
						
						%Rectangular pulse
						%phi_matrix(input_indexes,duration_stimulus)=phi_matrix(input_indexes,duration_stimulus)+obj.stimAmplitude;
						
						if pulses==0 %Control of markers and delivery of two pulses
							pulses=1;
							marker=1;
						else
							pulses=3;
							marker=3;
							flag_pause=1;
							obj.timePauseStart=n*obj.dt;
							obj.COUNT_STIMULI=obj.COUNT_STIMULI+1;
							%fprintf("start pause \n")
						end
						%Add the stim counter
						flags(2)=flags(2)+1;
					end
                        
				else
					%SHAM
					if pulses==0
						pulses=1;
						marker=2;
					else
						pulses=3;
						marker=4;
					end
				end
            end
            
            %Control for number of pulses (post-pulse procedures)
            if (flag_hold==1 && obj.phase>pi)
                if pulses==1
					%ready for the scond consecutive pulse
                    if (filtro_delta<obj.SOth)
                        flag_search=1;
                        pulses=2;
                    end
                elseif pulses==3
					%ready to find an SO again
                    flag_hold=0;
                    pulses=0;
                end
            end
            if (n*obj.dt-obj.timePauseStart)>(obj.pauseTime/obj.dt) && flag_pause
                flag_pause=0;
                flag_hold=0;
                pulses=0;
                fprintf('End pause \n')
            end
            if isempty(marker)
				marker=0;
			end
			stimulationFlags=[flags(1:9),flag_shot,flag_pause,flag_hold,flag_search,pulses];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %Phase stimulation
         function [obj,shape_now,marker,stimulationFlags]=phaseLock(obj,flags,n)
			currentSample=flags(1);
			old_marker=flags(7);
			V=flags(9);
			flag_shot=flags(10);
			flag_pause=flags(11);
			flag_hold=flags(12);
			flag_search=flags(13);
			pulses=flags(14);
			EPSILON_ANGLE=0.005*pi;
			marker=[];
			stimulus=zeros(obj.Nx*obj.Ny,1);
			STIM_SHAM=1;
            obj.previous_delta(2:end)=obj.previous_delta(1:end-1);
            obj.previous_delta(1)=V-obj.meanPhi;
            shape_now=0;
            nsamples=5000; %0.4 second para h=1e-4
            %filtering
            [filtro_delta,obj.previous_delta_filter]=util.filtrar(obj.numDelta,obj.denDelta,obj.previous_delta,obj.previous_delta_filter); %filtrar
            obj.x_now(2:3)=obj.x_now(1:2);
            obj.y_now(2:3)=obj.y_now(1:2);
            %Rectifiying
            if filtro_delta>=0
                obj.x_now(1)=filtro_delta;
            else
                obj.x_now(1)=-filtro_delta;
            end
			%New envelope detection method
			if obj.envelope_flag==0 && obj.envelope_count>0 && n>nsamples
				%Calculate a new tolerance afte nsamples of simulated data
				obj.envelope_tol=obj.envelope_accum/obj.envelope_count;
				obj.envelope_flag=1;
                fprintf('Tolerance envelope %.4e\n',obj.envelope_tol);
			end
			if obj.x_now(1)<obj.envelope_tol
				%Set the search for peak value to 0.01 if the signal is below the tolerance
				obj.peak_value=obj.ENVELOPE_OFFSET;
			elseif abs(obj.x_now(1)-obj.x_now(2))<=obj.envelope_tol && obj.x_now(1)>=obj.x_now(2)
				%when the signal arrives to a peak
				if obj.x_now(3)<obj.peak_value && obj.x_now(1)>obj.peak_value 
					%If the point is higher than the previous value, search for the peak value 
					obj.envelope_wait=true;
					if n<nsamples
						%Add the derivative to calculate the average tolerance value for zero derivative
						obj.envelope_accum=obj.envelope_accum+abs(obj.x_now(1)-obj.x_now(2));
						obj.envelope_count=obj.envelope_count+1;
					end
				else
					%If the point is not the peak, save the value to search for higher values than it
					obj.peak_value=obj.x_now(1);
				end
            elseif obj.envelope_wait==true && obj.x_now(1)<obj.x_now(2)
                %If the value is lesser than the previous one, the peak was reached, then save that value as the envelope 
                obj.previous_envelope=obj.x_now(2);
                obj.envelope_wait=false;
			end

			envolvente=obj.previous_envelope;
			if envolvente<obj.ENVELOPE_OFFSET %envelope lesser than 0.01
				envolvente=obj.ENVELOPE_OFFSET; %envelope set to 0.01
				%avoid the stimuli delivery until next oscillation			
			end
            if envolvente<2*obj.ENVELOPE_OFFSET
                flag_hold=1;
            end
			%Save the actual envelope
			obj.previous_envelope=envolvente;
			%delta(i)=filtro_delta;
			%Normalización
			eeg_normalizado=filtro_delta/envolvente;
			
			
			
			%Shifting phase -270 degrees (or forwarding 90 degrees)
			obj.shiftphase_eeg(2:end)=obj.shiftphase_eeg(1:end-1);
			obj.shiftphase_eeg(1)=eeg_normalizado;
			[shiftphase,obj.previous_shift_filter]=util.filtrar(obj.numShift,obj.denShift,obj.shiftphase_eeg,obj.previous_shift_filter);
			
			if shiftphase>=0
				obj.y_now(1)=shiftphase;
			else
				obj.y_now(1)=-shiftphase;
			end
			

			%New envelope detection method
			if obj.envelope_flag_shift==0 && obj.envelope_count_shift>0 && n>nsamples
				obj.envelope_tol_shift=obj.envelope_accum_shift/obj.envelope_count_shift;
				obj.envelope_flag_shift=1;
                fprintf('Tolerance envelope shift %.4e\n',obj.envelope_tol_shift);
			end
			if obj.y_now(1)<obj.envelope_tol_shift
				obj.peak_value_shift=obj.ENVELOPE_OFFSET;
			elseif abs(obj.y_now(1)-obj.y_now(2))<=obj.envelope_tol_shift && (obj.y_now(1)-obj.y_now(2))>=0
				if obj.y_now(3)<obj.peak_value_shift && obj.y_now(1)>obj.peak_value_shift
				    obj.envelope_wait_shift=true;	
					if n<nsamples
						obj.envelope_accum_shift=obj.envelope_accum_shift+abs(obj.y_now(1)-obj.y_now(2));
						obj.envelope_count_shift=obj.envelope_count_shift+1;
					end
				else
					obj.peak_value_shift=obj.y_now(1);
				end
            elseif obj.envelope_wait_shift==true && obj.y_now(1)<obj.y_now(2)
                obj.previous_envelope_shift=obj.y_now(2);
                obj.envelope_wait_shift=false;
                
			end

			envolvente_shift=obj.previous_envelope_shift;
			if envolvente_shift<0.1%The envelope of the shifted phase signal at leat must be 0.1
				envolvente_shift=0.1; %Set the envelope to 0.1
			end
			%Save the actual envelope
			obj.previous_envelope_shift=envolvente_shift;

			eeg_normalizado_90=shiftphase/envolvente_shift;
              
			%Calculating the current phase with the arctang
			%the shifted-phase signal is a sine
			%the original signal is takin in account as the cosine
			fase=atan(eeg_normalizado_90/(eeg_normalizado+1e-12));
			
			faseraw=fase; %save for comparisons
			
			%Translate from [-pi/2, pi/2] domain to [0, 2pi] range
			if obj.shiftphase_eeg(1)>0 && obj.shiftphase_eeg(2)<0 %cosine positively crossing zero 
				%aumento_fase=pi/2; python
				obj.phaseOffset=pi/2;
			elseif obj.shiftphase_eeg(2)>0 && obj.shiftphase_eeg(1)<0 %cosine neagtively crossing zero
				%aumento_fase=3/2*pi; python
				obj.phaseOffset=pi*3/2;
			end
			%sumar el desfase correspondiente segun el
			%cruce por cerostart
			fase=fase+obj.phaseOffset;
			delta=filtro_delta;
			%Si la fase resulta mayor a 2pi mantener en esa cota
			if fase>2*pi
				obj.phase=2*pi;
			elseif fase<0
				obj.phase=0;
			else
				obj.phase=fase;
			end
			
			%fprintf('%.6e\t %.6e\t %.6e\t %.6e \n',filtro_delta,envolvente,eeg_normalizado,obj.phase);  
			if flag_hold==0
				flag_search=1;
				%fprintf('flag_search set \n');
			end
			
			if (flag_search==1 && obj.phase>=obj.targetPhase-EPSILON_ANGLE && obj.phase<=obj.targetPhase+EPSILON_ANGLE)		        
				if flag_hold==0 && flag_shot==0
                    flag_shot=1;
                    flag_hold=1;  
                    
                    flags(1)=n;
                    flags(3)=1;
                            
                    %Rectangular pulse
                    %phi_matrix(input_indexes,duration_stimulus)=phi_matrix(input_indexes,duration_stimulus)+obj.stimAmplitude;
                    marker=1;
                    %Add the stim counter
                    flags(2)=flags(2)+1;
                    flag_search=0;  %Deactivate for not search the angle again until the first pulse is delivered
                end
			end

			
            %Control for search in the next low-frequency oscillation
            if (flag_hold==1 && obj.phase>1.9*pi)
                %fprintf('Reset at time: %.2f\n',n*1e-4);
                flag_hold=0;
                flag_shot=0;
            end
            if isempty(marker)
				marker=0;
			end
			stimulationFlags=[flags(1:9),flag_shot,flag_pause,flag_hold,flag_search,pulses];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        function setSOth(obj,SOth)
            %Update SO threshold
            obj.SOth=SOth;
        end
        
    end
end


