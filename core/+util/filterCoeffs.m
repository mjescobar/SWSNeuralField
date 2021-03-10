function [num,den,num1,den1]=filterCoeffs(freq,fs)
	%Conjugate poles near the unit circle at the desired frequency
    den=conv([1 -0.9999*(cos(freq*2*pi/fs)+1j*sin(freq*2*pi/fs))],[1 -0.9999*(cos(freq*2*pi/fs)-1j*sin(freq*2*pi/fs))]);
    %Numerator 
    num=[1, 1, -1, -1]*1e-4;
    
    den1=conv(den,[1 -0.999]);
    %Conjugate zeros with radius=1/2 very near to angle pi/2 
    numaux=conv([1 0.5*(cos(0.99*pi/2)+1j*sin(0.99*pi/2))],[1 0.5*(cos(0.99*pi/2)-1j*sin(0.99*pi/2))]);
    %numaux1=conv([1 0.9*(cos(freq*2*pi/fs)+1j*sin(freq*2*pi/fs))],[1 0.9*(cos(freq*2*pi/fs)-1j*sin(freq*2*pi/fs))]);
    %num1=conv(numaux,numaux1);
    num1=conv(numaux,[1 0.9]);
    %Add the pole at zero for the invert-notch filter
    den=[den,0];
    
    factor=sum(den1)/sum(num1);
    num1=num1*factor/2;
    
end
