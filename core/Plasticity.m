classdef Plasticity
	properties
		HGamma
		H
		ratesA
		ratesB
		omega
		Phifft
		N
		nfft
		tp
		gamma
		range
	end
	methods
		function obj=Plasticity(config)
			%Values and calculations from 'Assadzadeh, 2018'
			%obj.ratesA=[-1.01;-1.01;-1.25;-1;-0.4;-1.25;-1.25;-1.25];
			obj.ratesA=config.ratesA;
			%obj.ratesB=[10,1.6,0.8,2.7,0.7,1.1,0.55,0.45];
			obj.ratesB=config.ratesB;
			obj.tp=config.tp;
			obj.omega=config.omega;
			obj.N=config.Nx*config.Ny;
			obj.Phifft=zeros(5,config.samplingFrequency*2);
			obj.nfft=config.samplingFrequency;
			obj.gamma=config.gamma;
			obj.range=config.range;
			fprintf('Plasticity: Instatiation completed\n')
		end

		function obj=init(obj)
			%Calculate transfer function spectrums
			obj.H=obj.calculateH();%Verify this
			Gamma=obj.calculateGamma(obj.gamma,obj.range);
			obj.HGamma=conj(obj.H(1,:)).*conj(Gamma);
		end
		
		function H=calculateH(obj)
			%STDP spectrum
			H=obj.ratesB'.*(obj.ratesA*obj.tp*(1+1j*obj.omega)+obj.tp*(1-1j*obj.omega))./(1+(obj.omega*obj.tp).^2);
		end
		
		function [mfft,means]=calculatefft(obj,Phi)
			%Populations spectrums
			mfft=zeros(5,obj.nfft);
			means=zeros(5,1);
			for m=0:4
				[Phi_pop,meanPhi]=obj.zScore(Phi(m*obj.N+1:(m+1)*obj.N,:));
				means(m+1)=meanPhi;
				sfft=fft(Phi_pop);
				mfft(m+1,:)=sfft(1:obj.nfft);
			end
			
		end
			
		function [z,meanx]=zScore(obj,x)
			%z-score of 
			meanx=mean(x(:));
			z=(mean(x,1)-meanx)./std(x(:));
		end
		
		function Gamma=calculateGamma(obj,gamma,range)
			k=(2*pi*3e9)./obj.omega;
			Gamma=((1-1j*obj.omega./gamma).^2+k.^2*range^2).^(-1);
		end
		
		function integ=integrand(obj,phia,phib,indexH,flagGamma)
			if flagGamma==1
				integ=real(phia.*conj(phib).*obj.HGamma);
			else
				integ=real(phia.*conj(phib).*conj(obj.H(indexH,:)));
			end
			
		end
		
		function integrall=integrall(obj,integrand)
			integrall=integrand(1)*(obj.omega(2)-obj.omega(1))/(2*pi);
			for i=2:length(integrand)
				integrall=integrall+integrand(i)*(obj.omega(i)-obj.omega(i-1))/(2*pi);
			end
		end
		
		function [strengthRatios,means]=update(obj,Phi)
			[mfft,means]=obj.calculatefft(Phi);
			strengthRatios=zeros(8,1);
			%s_ee
			strengthRatios(1)=obj.integrall(obj.integrand(mfft(1,:),mfft(1,:),1,0));
			%s_ei
			strengthRatios(2)=obj.integrall(obj.integrand(mfft(1,:),mfft(2,:),2,0));
			%s_es
			strengthRatios(3)=obj.integrall(obj.integrand(mfft(1,:),mfft(4,:),3,0));
			%s_re
			strengthRatios(4)=obj.integrall(obj.integrand(mfft(3,:),mfft(1,:),4,0));
			%s_rs
			strengthRatios(5)=obj.integrall(obj.integrand(mfft(3,:),mfft(4,:),5,0));
			%s_se
			strengthRatios(6)=obj.integrall(obj.integrand(mfft(4,:),mfft(1,:),6,0));
			%s_sr
			strengthRatios(7)=obj.integrall(obj.integrand(mfft(4,:),mfft(3,:),7,0));
			%s_sn
			strengthRatios(8)=obj.integrall(obj.integrand(mfft(4,:),mfft(5,:),8,0));
		end
	end
end
