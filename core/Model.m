classdef Model
    %"""Solver class, to define paramaeters of the model and to tag results"""
    properties
        connectionMatrix %Connections matrix
        connectionMatrixDelay %Connections with delay
        connectionMatrixNodelay %Connections without delay
        plasticityMatrix %Plasticity changing strength ratios
        plasticityVector %changes in plasticity.
        connectionVector %changes in plasticity.
        delayMatrixMask %maskMatrix
        Nx %"Horizontal" size of a population 
        Ny %"Vertical"size of a population
        L %Total nodes per population
        W %Total connection matrix size M*Lx*Ly
        %y=Ax+Bu
        A %Matrix A 
        B %Matrix B
        alpha
        beta
        gamma %damping factor
        range %axonal range
        t0 %Corticothalamic delay
        Qmax %Max firing rate
        theta %Soma potential threshold
        sigma_rho %slope sigmoid->deviation of the normal distribution of the firing response
        inverseSigma_rho
        stencilAAS
        ee_idx
        ei_idx
        es_idx
        re_idx
        rs_idx
        se_idx
        sr_idx
        sn_idx
    end
    
    methods
		%Constructor
        function obj=Model(config)
            obj.connectionMatrix=config.connectionMatrix;
            obj.connectionMatrixDelay=config.connectionMatrixDelay;
			obj.connectionMatrixNodelay=config.connectionMatrixNodelay;
            obj.plasticityVector=config.initialChangeRates;
            obj.connectionVector=config.initialStrengths;
            obj.Nx=config.Nx;
            obj.Ny=config.Ny;
            obj.L=obj.Nx*obj.Ny;
            obj.W=5*obj.L;
            obj.alpha=config.alpha;
            obj.beta=config.beta;
            %%%Old way to insert the time delay in the connections
            delay=[0;0;1;1;0;1;0;0];
            delay_connections=util.eirsConnections(delay);
            obj.delayMatrixMask=util.expandConnections(obj.Nx,obj.Ny,delay_connections);
            [obj.ee_idx,obj.ei_idx,obj.es_idx,obj.re_idx,obj.rs_idx,obj.se_idx,obj.sr_idx,obj.sn_idx]=util.indicesConnections(obj.Nx,obj.Ny);
			%obj.connectionMatrixDelay=obj.connectionMatrix.*obj.delayMatrix;
			%obj.connectionMatrixNodelay=obj.connectionMatrix.*(1-obj.delayMatrix);
			%EIRS MODEL
            obj.A=[0, 1; -obj.alpha*obj.beta, -obj.alpha-obj.beta]';
            obj.B=[0; obj.alpha*obj.beta]';
            obj.t0=config.t0; %seconds
            obj.gamma=config.gamma;
            obj.range=config.range;
            obj.Qmax=config.Qmax;
            obj.theta=config.theta;
            obj.sigma_rho=config.sigma_rho;
            obj.inverseSigma_rho=1/obj.sigma_rho;
            stencilAAS=zeros(obj.W,1);
            if config.AASconnections==1
				stencilAAS(4*obj.L+1:end,1)=1;
			elseif config.AASconnections==2
				stencilAAS(3*obj.L+1:4*obj.L,1)=1;
			end
            obj.stencilAAS=stencilAAS;
            
            fprintf('Model: Instantiatiation completed \n')
            
        end
        
        %% Update connections strengths by plasticity
        function obj=updateConnections(obj,twindow,tau,strengthsRates)
			obj.plasticityVector=strengthsRates;
			obj.connectionVector=obj.connectionVector+twindow/tau.*strengthsRates;
			[obj.connectionMatrix,obj.connectionMatrixDelay, obj.connectionMatrixNodelay]=util.connectionsUpdate(obj.connectionVector,...
			obj.connectionMatrix, obj.connectionMatrixDelay, obj.connectionMatrixNodelay,obj.ee_idx,obj.ei_idx,obj.es_idx,obj.re_idx,obj.rs_idx,obj.se_idx,obj.sr_idx,obj.sn_idx);
        end
        
        %% Sigmoid
        function Q=firingResponse(obj,V)
			%Q=obj.Qmax./(1+exp(-(V-obj.theta)*obj.inverseSigma_rho));
			Q=firingResponse(V,obj.Qmax,obj.theta,obj.inverseSigma_rho);
		end
		
		%% Sigmoid approx
		function Q=firingResponseApprx(obj,V)
			Q=obj.Qmax*(1+util.expTruncSerie(-(V-obj.theta)*obj.inverseSigma_rho)).^(-1);
        end
        %% Linear activation
         %function firing_response=linear(obj,V,Qmax,Q0,a)
            %"""Sigmoid function for firing response"""
         %   firing_response=Q0+a.*V;
          %  if firing_response>Qmax
          %      firing_response=Qmax;
          %  end
        %end
		%% Dynamic matrices
		function obj=updateDerivateMatrices(obj)
			obj.A=[0, 1; -obj.alpha*obj.beta, -obj.alpha-obj.beta]';
            obj.B=[0; obj.alpha*obj.beta]';
        end
        %% Change values of alpha and beta
        function obj=setalphabeta(obj,alpha,beta)
			 obj.alpha=alpha;
			 obj.beta=beta;
         end
         
          %% Set connections strengths by force
        function obj=setConnections(obj,newStrengths)
			newConnections=util.eirsConnections(newStrengths);
            obj.connectionMatrix=util.expandConnections(obj.Nx,obj.Ny,newConnections);
			obj.connectionMatrixDelay=obj.connectionMatrix.*obj.delayMatrix;
			obj.connectionMatrixNodelay=obj.connectionMatrix.*(1-obj.delayMatrix);
		end
    end
end
