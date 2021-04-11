function [output_matrix,matrix_delay,matrix_nodelay]=expand_connections(Nx,Ny,strengths_matrix)
	%"""Expands each strenght between populations as an intra-population identity matrix """
	size_M=size(strengths_matrix);
	L=Nx*Ny;
	W=size_M(1)*L;
	identity=eye(L);
	connections_matrix=zeros(W);
	matrix_delay=zeros(W);
	matrix_nodelay=zeros(W);
	for i =0:size_M(1)-1
		for j=0:size_M(2)-1
			connections_matrix(i*L+1:(i+1)*L,j*L+1:(j+1)*L)=identity*strengths_matrix(i+1,j+1);
			if any(i==[0,1]) && j==3
				matrix_delay(i*L+1:(i+1)*L,j*L+1:(j+1)*L)=identity*strengths_matrix(i+1,j+1);
			elseif any(i==[2,3]) && j==0
				matrix_delay(i*L+1:(i+1)*L,j*L+1:(j+1)*L)=identity*strengths_matrix(i+1,j+1);
			else
				matrix_nodelay(i*L+1:(i+1)*L,j*L+1:(j+1)*L)=identity*strengths_matrix(i+1,j+1);
		end
	end
	output_matrix=connections_matrix;
end
