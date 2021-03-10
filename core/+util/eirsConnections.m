 function strengths_matrix=eirsConnections(eirs_strengths)
		%"""Construct a matrix with the relevant strenghts"""
		%%% EIRS weights (in some papers the gains are indicated Gab=rho*vab, and rho is the slope of the sigmoid ~340)
		% [vee, vei, ves, vie, vii, vis, vir, vre, vrs, vse, vsr, vsn]
		% vee==vie
		% vei==vii
		% ves==vis
		% [vee, vei, ves, vre, vrs, vse, vsr, vsn]
		strengths_matrix=zeros(5,5);
		strengths_matrix(1,1)=eirs_strengths(1);
		strengths_matrix(1,2)=eirs_strengths(2);
		strengths_matrix(1,4)=eirs_strengths(3);
		strengths_matrix(2,1)=eirs_strengths(1);
		strengths_matrix(2,2)=eirs_strengths(2);
		strengths_matrix(2,4)=eirs_strengths(3);
		strengths_matrix(3,1)=eirs_strengths(4);
		strengths_matrix(3,4)=eirs_strengths(5);
		strengths_matrix(4,1)=eirs_strengths(6);
		strengths_matrix(4,3)=eirs_strengths(7);
		strengths_matrix(4,5)=eirs_strengths(8);
end
