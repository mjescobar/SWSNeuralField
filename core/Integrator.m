classdef Integrator

	properties
		A
		Atilde
		B
		Btilde
		Atilde_aas
		Btilde_aas
		h
		propagator
	end

	methods
		function obj=Integrator_old(config,model)
			obj.h=config.h;
			obj.A=model.A;
			obj.B=model.B;
			obj.Atilde=(eye(2)+obj.h*obj.A);
			obj.Btilde=(obj.h*obj.B);
			obj.propagator=Propagator(config,model);
		end
		function obj=Integrator(config,model,modelAAS)
			obj.h=config.h;
			obj.A=model.A;
			obj.B=model.B;
			obj.Atilde=(eye(2)+obj.h*obj.A);
			obj.Btilde=(obj.h*obj.B);
			obj.Atilde_aas=(eye(3)+obj.h*modelAAS.A);
			obj.Btilde_aas=obj.h*modelAAS.B;
			obj.propagator=Propagator(config,model);
		end
		
		function xnplusone=stepEuler(obj,x,u)
			%Euler integration method
			%xnplusone=x+obj.h*(x*obj.A+u*obj.B);
			%xnplusone=x*obj.Atilde+u*obj.Btilde;
			xnplusone=prodsum(obj.Atilde,x,obj.Btilde,u);
		end
		function [xnplusone,xnplusone_aas]=stepEulerAAS(obj,x,u,x_aas,u_aas)
			%Euler integration method
			%xnplusone=x+obj.h*(x*obj.A+u*obj.B);
			%xnplusone=x*obj.Atilde+u*obj.Btilde;
			xnplusone=prodsum(obj.Atilde,x,obj.Btilde,u);
			%xnplusone_aas=x_aas*obj.Atilde_aas+u_aas*obj.Btilde_aas;
			xnplusone_aas=prodsum(obj.Atilde_aas,x_aas,obj.Btilde_aas,u_aas);
		end
		function xnplusone=stepHeun(obj,x,u)
			%Heun integration method
			dx=(x*obj.A+u*obj.B);
			xtilde=x+obj.h*dx;
			xnplusone=x+obj.h/2*(dx+xtilde*obj.A+u*obj.B);
		end	
	end
end
