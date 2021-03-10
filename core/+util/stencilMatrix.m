function S=stencilMatrix(Nx,Ny,p2)
W=Nx*Ny;
S=zeros(W);
nonZeroIndexes=zeros(W,4);
for m=1:Nx %filas
    for n=1:Ny %columnas
        %Left
        node=(m-1)*Ny+n;
        if n==1
			nonZeroIndexes(node,1)=m*Ny;
        else
           nonZeroIndexes(node,1)=(m-1)*Ny+n-1;
        end
        
         %Up
        if m==1
			nonZeroIndexes(node,2)=(Nx-1)*Ny+n; 
        else
			nonZeroIndexes(node,2)=(m-2)*Ny+n;
        end
        
        %Right
        if n==Ny
			nonZeroIndexes(node,3)=(m-1)*Ny+1;  
        else
			nonZeroIndexes(node,3)=(m-1)*Ny+n+1;
        end
        
        %Down
        if m==Nx
			nonZeroIndexes(node,4)=n;   
        else
			nonZeroIndexes(node,4)=m*Ny+n;
        end
        
    end
end
for j=1:W
    for k=1:4
    S(j,nonZeroIndexes(j,k))=p2;
    end
end
