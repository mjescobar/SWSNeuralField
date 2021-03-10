function [ee,ei,es,re,rs,se,sr,sn]=indicesConnections(Nx,Ny)
	%Return the indices of the connection matrix for each strength value
	M=5; %populations
	L=Nx*Ny; %Extent of the matrix of each population
	W=M*L; %Extent of the entire connection matrix
	connections_matrix=zeros(W);
	identity=eye(L);
	n=1;
	strengths_matrix=zeros(5,5);
	strengths_matrix(1,1)=1;
	strengths_matrix(1,2)=2;
	strengths_matrix(1,4)=3;
	strengths_matrix(2,1)=1;
	strengths_matrix(2,2)=2;
	strengths_matrix(2,4)=3;
	strengths_matrix(3,1)=4;
	strengths_matrix(3,4)=5;
	strengths_matrix(4,1)=6;
	strengths_matrix(4,3)=7;
	strengths_matrix(4,5)=8;
	for i =0:M-1
		for j=0:M-1
			connections_matrix(i*L+1:(i+1)*L,j*L+1:(j+1)*L)=identity*strengths_matrix(i+1,j+1);
		end
	end
	ee=find(connections_matrix==1);
	ei=find(connections_matrix==2);
	es=find(connections_matrix==3);
	re=find(connections_matrix==4);
	rs=find(connections_matrix==5);
	se=find(connections_matrix==6);
	sr=find(connections_matrix==7);
	sn=find(connections_matrix==8);
end
