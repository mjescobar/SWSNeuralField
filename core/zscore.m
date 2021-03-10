function zx=zscore(x)
	zx=(mean(x,1)-mean(x))./std(x);
end
