function filename=buildFilename(simulationPoint,stimShape,stimFrequency, stimAmplitude, stimPulseDuration, stimAmplitudeSD, stimFrequencySD, plasticity_onoff)


	if stimAmplitudeSD==0
		sdA='OFF';
	else
		sdA='ON';
	end
	if stimFrequencySD==0
		sdF='OFF';
	elseif stimFrequencySD==100
		sdF='Poisson';
	else
		sdF='ON';
	end
	%shape: 0: rectangular, 1: rise ramp, 2: decrease ramp, 3: triangle,
     %4: gaussian, 5: Hamming, 6:sinc2, 7:sin, 8:cos
	if stimShape==0
		shape='rectangular';
	elseif stimShape==1
	shape='riseRamp';
	elseif stimShape==2
	shape='decreaseRamp';
	elseif stimShape==3
	shape='triangle';
	elseif stimShape==4
	shape='Gaussian';
	elseif stimShape==5
	shape='hamming';
	elseif stimShape==6
	shape='sinc';
	elseif stimShape==7
	shape='sin';
	elseif stimShape==8
	shape='sin2';
	elseif stimShape==9
	shape='cos';
	elseif stimShape==10
	shape='expLongTail';
	elseif stimShape==11
	shape='expShortTail';
	elseif stimShape==12
	shape='trapece';
	elseif stimShape==13
	shape='halfTrapece';
	elseif stimShape==14
	shape='exponential';
	end
	if stimAmplitude==0
	filename=sprintf('%s-%s',simulationPoint,'baseline');
	else
		if plasticity_onoff==0
		filename=sprintf('%s-%s-F%.2f-A%.2e-D%.2e-sdA-%s-sdF-%s',simulationPoint,shape,stimFrequency,stimAmplitude,stimPulseDuration,sdA,sdF);
		else
		filename=sprintf('%s-p-%s-F%.2f-A%.2e-D%.2e-sdA-%s-sdF-%s',simulationPoint,shape,stimFrequency,stimAmplitude,stimPulseDuration,sdA,sdF);
		end
	end

end
