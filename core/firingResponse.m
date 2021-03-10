function Q=firingResponse(V,Qmax,theta,inv_sigma)
	Q=Qmax./(1+exp(-(V-theta)*inv_sigma));
end
