function [y,y1]=filtrar(b,a,x,y1)
y=b*x-a(2:end)*y1;
y1(2:end)=y1(1:end-1);
y1(1)=y;
end
