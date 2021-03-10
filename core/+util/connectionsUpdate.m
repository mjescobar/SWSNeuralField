function [connection_matrix,matrix_delay,matrix_nodelay]=connectionsUpdate(strengths_vector,connection_matrix,matrix_delay,matrix_nodelay,ee,ei,es,re,rs,se,sr,sn)
	connection_matrix(ee)=strengths_vector(1);
	connection_matrix(ee)=strengths_vector(2);
	connection_matrix(ee)=strengths_vector(3);
	connection_matrix(ee)=strengths_vector(4);
	connection_matrix(ee)=strengths_vector(5);
	connection_matrix(ee)=strengths_vector(6);
	connection_matrix(ee)=strengths_vector(7);
	connection_matrix(ee)=strengths_vector(8);
	%Delay
	matrix_delay(es)=strengths_vector(3);
	matrix_delay(re)=strengths_vector(4);
	matrix_delay(se)=strengths_vector(6);
	%No delay
	matrix_nodelay(ee)=strengths_vector(1);
	matrix_nodelay(ei)=strengths_vector(2);
	matrix_nodelay(rs)=strengths_vector(5);
	matrix_nodelay(sr)=strengths_vector(7);
	matrix_nodelay(sn)=strengths_vector(8);
end
