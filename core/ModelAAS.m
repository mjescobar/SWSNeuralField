%AAS model
classdef ModelAAS
	properties
		nu_vm
		nu_mv
		nu_vc
		nu_vh
		tau_m
		tau_v
		xi
		mu
		c0
		Cphase
		Qm
		Qv
		Vaas
		A
		B
		Qmax_AAS
		inverseSigma_rho
		theta	
	end
	
	methods
		function obj=ModelAAS(config)
			%Set constants
			obj.nu_vm=config.nu_vm;
			obj.nu_mv=config.nu_mv;
			obj.nu_vc=config.nu_vc;
			obj.nu_vh=config.nu_vh;
			obj.xi=config.xi;
			obj.c0=config.c0;
			obj.mu=config.mu;
			obj.tau_m=config.tau_m;
			obj.tau_v=config.tau_v;
			obj.Cphase=config.Cphase;
			obj.Qmax_AAS=100;
			obj.theta=config.theta;
			obj.inverseSigma_rho=1/config.sigma_rho;
			%AAS Model
            obj.A=[-1/config.tau_v,0,0;0,-1/config.tau_m,0;0,0,-1/config.xi]';
            obj.B=[1,0,0;0,1,0;0,0,1]';
        end
			
		function Q=firingResponseAAS(obj,V)
			%Q=obj.Qmax_AAS./(1+exp(-(V-obj.theta)*obj.inverseSigma_rho));
			Q=firingResponse(V,obj.Qmax_AAS,obj.theta,obj.inverseSigma_rho);
		end
		

		function D=sleepDrive(obj,time,H)
			D=obj.nu_vc*(obj.c0+cos(2*pi*time/(24*3600)+obj.Cphase))+obj.nu_vh*H;
		end
		
		function u=vector_u(obj,Q_aas,D,Ach)
			u=[(obj.nu_vm*Q_aas(2)+D)/obj.tau_m,(obj.nu_vm*Q_aas(1)+Ach)/obj.tau_m,obj.mu*Q_aas(2)/obj.xi];
			
		end
	end
end

