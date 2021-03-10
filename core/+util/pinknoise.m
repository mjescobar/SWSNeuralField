function [x,randreg]=pinknoise(randreg)
%randreg=zeros(5,1)
    %Author: http://www.ridgerat-tech.us/tech/pinkalg.htm
    Ngen = 3;
    am=sqrt(2)*217; %(217~1/E{randreg})
    %E{randreg}=(pv(2)-pv(1))*(av(1))+(pv(3)-pv(2))*(av(1)+av(2))+(1-pv(3))*(av(1)+av(2)+av(3));
    av = [ 4.6306e-003  5.9961e-003  8.3586e-003 ];
    pv = [ 3.1878e-001  7.7686e-001  9.7785e-001  ];
    if sum(randreg)==0
        for igen=1:Ngen
          randreg(igen)=av(igen)*am*(rand(1,1)-0.5);
        end
    end
     rv = rand(1,1);

      % Update each generator state per probability schedule
      for igen=1:Ngen
        if  rv > pv(igen)
          randreg(igen) = av(igen)*am*(rand(1,1)-0.5);
        end
      end

      % Signal is the sum of the generators
      x=sum(randreg);