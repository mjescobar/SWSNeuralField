classdef Monitor
	properties
	filename
	indexesV
	indexesQ
	indexesPhi
    %indexesStim
	fileQ
	fileV
	filePhi
    fileStim
    filePhase
    fileStrengths
    %batch storage
    nrows
	end
	methods
		%Constructor
		function obj=Monitor(config,filename,signals)
			obj.filename=filename;
			obj.indexesV=[];
			obj.indexesQ=[];
			obj.indexesPhi=[];
			obj.nrows=config.storageRows;
			fprintf('Monitor: Starting new monitor\n')
            %obj.indexesStim=1:config.Nx*config.Ny;
            %Select the indexes of the required nodes
			obj.fileStim=fopen([filename,'-Stim.txt'],'w+');
			for m=1:length(signals)
				str=split(signals{m},':');
				if str{3}=="all"
					indexes=1:config.Nx*config.Ny;
					if str(2)=="E"
						indexes=indexes;
					elseif str{2}=="I"
						indexes=indexes+config.Nx*config.Ny;
					elseif str{2}=="R"
						indexes=indexes+2*config.Nx*config.Ny;
					elseif str{2}=="S"
						indexes=indexes+3*config.Nx*config.Ny;
					elseif str{2}=="N"
						indexes=indexes+4*config.Nx*config.Ny;
					end
				else
					[indexes,n,errmsg]=sscanf(str{3},'%d');
					if min(indexes)<1 || max(indexes)>config.Nx*config.Ny
						fprintf('Error, indexes could not be negative or superior to %d',config.Nx*config.Ny);
					else			
						if str(2)=="E"
							indexes=indexes;
						elseif str{2}=="I"
							indexes=indexes+config.Nx*config.Ny;
						elseif str{2}=="R"
							indexes=indexes+2*config.Nx*config.Ny;
						elseif str{2}=="S"
							indexes=indexes+3*config.Nx*config.Ny;
						elseif str{2}=="N"
							indexes=indexes+4*config.Nx*config.Ny;
						end
					end
                end
                %open .txt files
                
				if str{1}=="V"
					obj.indexesV=[obj.indexesV,indexes];
					obj.fileV=fopen([filename,'-V.txt'], 'w+');
				elseif str{1}=="Q"
					obj.indexesQ=[obj.indexesQ,indexes];
					obj.fileQ=fopen([filename,'-Q.txt'], 'w+');
				elseif str{1}=="Phi"
					obj.indexesPhi=[obj.indexesPhi,indexes];
					obj.filePhi=fopen([filename,'-Phi.txt'], 'w+');
				else
					fprintf('Monitor: missunderstood expression in signals %d \n',m);
				end
			end
			
			fprintf('Monitor: Instatiation completed for %s\n',filename)
		end
		
		%Save data in a text file
		function flag=save(obj,V,Q,phi,stim,marker,currentTime)
            fprintf(obj.fileStim,'%.9e\t',stim);
            fprintf(obj.fileStim,'%d\t%.4e',marker,currentTime);
            fprintf(obj.fileStim,'\n');
			if ~isempty(obj.indexesV)
			fprintf(obj.fileV,'%.9e \t',V(obj.indexesV));
			fprintf(obj.fileV,'\n');
			end
			if ~isempty(obj.indexesQ)
			fprintf(obj.fileQ,'%.9e \t',Q(obj.indexesQ));
			fprintf(obj.fileQ,'\n');
			end
			if ~isempty(obj.indexesPhi)
			fprintf(obj.filePhi,'%.9e \t',phi(obj.indexesPhi));
			fprintf(obj.filePhi,'\n');
			end
		end
		
		%Save batch data in a text file
		function flag=savebatch(obj,V,Q,phi,stim,marker,currentTime)
			for i=1:obj.nrows
				fprintf(obj.fileStim,'%.9e\t',stim(:,i));
				fprintf(obj.fileStim,'%d\t%.4e',marker(:,i),currentTime(:,i));
				fprintf(obj.fileStim,'\n');
			
				if ~isempty(obj.indexesV)
				fprintf(obj.fileV,'%.9e \t',V(obj.indexesV,i));
				fprintf(obj.fileV,'\n');
				end
				if ~isempty(obj.indexesQ)
				fprintf(obj.fileQ,'%.9e \t',Q(obj.indexesQ,i));
				fprintf(obj.fileQ,'\n');
				end
				if ~isempty(obj.indexesPhi)
				fprintf(obj.filePhi,'%.9e \t',phi(obj.indexesPhi,i));
				fprintf(obj.filePhi,'\n');
				end
			end
		end
		
		function flag=saveWithoutStim(obj,V,Q,phi)
			if ~isempty(obj.indexesV)
			fprintf(obj.fileV,'%.9e \t',V(obj.indexesV));
			fprintf(obj.fileV,'\n');
			end
			if ~isempty(obj.indexesQ)
			fprintf(obj.fileQ,'%.9e \t',Q(obj.indexesQ));
			fprintf(obj.fileQ,'\n');
			end
			if ~isempty(obj.indexesPhi)
			fprintf(obj.filePhi,'%.9e \t',phi(obj.indexesPhi));
			fprintf(obj.filePhi,'\n');
			end
		end
		%Close the txt files. (Call this function at the end of simulation)
		function close(obj)
			if ~isempty(obj.indexesV)
				fclose(obj.fileV);
			end
			if ~isempty(obj.indexesQ)
				fclose(obj.fileQ);
			end
			if ~isempty(obj.indexesPhi)
				fclose(obj.filePhi);
			end
			fclose(obj.fileStim);

		end
		
		%Phase
		function obj=createFilePhase(obj,filename)
			obj.filePhase=fopen([filename,'-Phase.txt'], 'w+');
			fprintf('Created file to store phase \n')
		end
		
		function obj=createFileStrengths(obj,filename)
			obj.fileStrengths=fopen([filename,'-Strengths.txt'], 'w+');
			fprintf('Created file to store strengths and mean values\n')
		end
		
		function savePhase(obj,Phase)
			if ~isempty(obj.filePhase)
				fprintf(obj.filePhase,'%.9f,%.9f,%.9f\n',Phase(1),Phase(2),Phase(3));
			end
		end
		
		function savePhasebatch(obj,Phase)
			if ~isempty(obj.filePhase)
				for i=1:obj.nrows
					fprintf(obj.filePhase,'%.9f\t',Phase(:,i));
					fprintf(obj.filePhase,'\n');
				end
			end
		end
		
		function saveStrengths(obj,strengths, meanPhi,ratios)
			if ~isempty(obj.fileStrengths)
					fprintf(obj.fileStrengths,'%.9f \t',strengths);
					fprintf(obj.fileStrengths,'%.9f \t',meanPhi);
					fprintf(obj.fileStrengths,'%.9f \t',ratios);
					fprintf(obj.fileStrengths,'\n');
			end
		end
		
		function closePhase(obj)
			if ~isempty(obj.filePhase)
				fclose(obj.filePhase);
			end
		end
		
		function closeStrengths(obj)
			if ~isempty(obj.fileStrengths)
				fclose(obj.fileStrengths);
			end
		end
	end
end
