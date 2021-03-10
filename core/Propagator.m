classdef Propagator
	properties
		dfact
		p2
		tenminus4p2
		twominus4p2
		expfactpos
		expfactneg
		coeficients1
		coeficients2
		I
		Sp2
	end
	methods
		%Constructor
		function obj=Propagator(config,model)
		dt2on12=config.h^2/12;
		obj.dfact=dt2on12*model.gamma^2;
		dt2ondx2=config.h^2/config.dx^2;
		p2=dt2ondx2*model.range^2*model.gamma^2;
		obj.tenminus4p2=10-4*p2;
		obj.twominus4p2=2-4*p2;
		obj.expfactneg=exp(-config.h*model.gamma);
        obj.expfactpos=exp(config.h*model.gamma);
		obj.I=eye(model.Nx*model.Ny);
		obj.Sp2=util.stencilMatrix(model.Nx,model.Ny,p2);
		obj.coeficients1=zeros(3,1);
		obj.coeficients2=zeros(2,1);
		obj.coeficients1(1,1)=obj.expfactpos;
		obj.coeficients1(2,1)=obj.tenminus4p2;
		obj.coeficients1(3,1)=obj.expfactneg;
		obj.coeficients2(1,1)=obj.twominus4p2;
		obj.coeficients2(2,1)=-obj.expfactneg;
		end
		%% Spatial propagation
			function phinplusone=propagate(obj,matrixQ3,matrixPhi2)
				%q_matrix_3 last 3 values
				%phi_matrix_3 last 2 values
				%time indexes
				%index 3: n-1, index 2: n, index 1: n+1
				%Indexes are obtained by cirshift([1,2,3],-1);
				%if deltat>gamma/2 && deltax>range/2
				%drive=obj.dfact*((obj.tenminus4p2*obj.I+obj.Sp2)*matrixQ3(:,2)+obj.expfactpos*obj.I*matrixQ3(:,1)+obj.expfactneg*obj.I*matrixQ3(:,3));
				drive=obj.dfact*(obj.Sp2*matrixQ3(:,2)+obj.I*matrixQ3*obj.coeficients1);
				%phinplusone=obj.expfactneg*((obj.twominus4p2*obj.I+obj.Sp2)*matrixPhi2(:,1)-obj.expfactneg*obj.I*matrixPhi2(:,2)+drive);
				phinplusone=obj.expfactneg*(obj.Sp2*matrixPhi2(:,1)+obj.I*matrixPhi2*obj.coeficients2+drive);
                %else
				%fprintf('gamma and range must be at maximum the double of delta t and delta x, respectively.')
				%end
			end
	end
end
