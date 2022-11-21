function params=subParameters(comm_params, sub_strains)

    totStrains=length(comm_params.Vs)/2;
    numStrains=length(sub_strains);
    
    params.seg_rate=comm_params.seg_rate;
    params.conj_rate=comm_params.conj_rate;
    params.S0=comm_params.S0;
    params.T=comm_params.T;
    params.d=comm_params.d;
    params.epsilon=comm_params.epsilon;
    
    params_species={};
    params_plasmids={};
    params_strains={};
    params_colors=zeros(2*numStrains,3);
    params_cs=zeros(1,2*numStrains);  %rhos
    params_Vs=zeros(1,2*numStrains);  %Vmax
    params_Ks=zeros(1,2*numStrains);  %Km
    
    params_MIC_AMC=zeros(1,2*numStrains);  
    params_MIC_MER=zeros(1,2*numStrains);  
    params_MIC_IMP=zeros(1,2*numStrains);  
    params_MIC_ERT=zeros(1,2*numStrains);  
    params_conj_freq=zeros(1,2*numStrains);  
    
    for icomm=1:numStrains
        
        this_strain=sub_strains{icomm};
       
       istrain_all = find(strcmp(comm_params.strains, this_strain));
       if ~isempty(istrain_all)
            istrain=min(istrain_all);
       
            %Plasmid-bearing (TC)
            params_cs(icomm)=comm_params.cs(istrain);
            params_Vs(icomm)=comm_params.Vs(istrain);
            params_Ks(icomm)=1;
            params_species{icomm}=comm_params.species{istrain};
            params_strains{icomm}=extractBefore(comm_params.strains{istrain},'_');
            params_plasmids{icomm}='TC';
            
            params_MIC_AMC(icomm)=comm_params.MIC_AMC(istrain);  
            params_MIC_MER(icomm)=comm_params.MIC_MER(istrain);  
            params_MIC_IMP(icomm)=comm_params.MIC_IMP(istrain);  
            params_MIC_ERT(icomm)=comm_params.MIC_ERT(istrain);  
            params_conj_freq(icomm)=comm_params.conj_freq(istrain);  
            
            if isfield(comm_params, 'colors')
                params_colors(icomm,:)=comm_params.colors(istrain,:); 
            else
                param_color(icomm,:)=[.5 .5 .5]; 
            end
            
            %Plasmid-bearing (TC)
            params_cs(icomm+numStrains)=comm_params.cs(istrain+totStrains);
            params_Vs(icomm+numStrains)=comm_params.Vs(istrain+totStrains);
            params_Ks(icomm+numStrains)=1;
            params_species{icomm+numStrains}=comm_params.species{istrain};
            params_strains{icomm+numStrains}=extractBefore(comm_params.strains{istrain},'_');
            params_plasmids{icomm+numStrains}='WT';
            
            
            params_MIC_AMC(icomm+numStrains)=comm_params.MIC_AMC(istrain+totStrains);  
            params_MIC_MER(icomm+numStrains)=comm_params.MIC_MER(istrain+totStrains);  
            params_MIC_IMP(icomm+numStrains)=comm_params.MIC_IMP(istrain+totStrains);  
            params_MIC_ERT(icomm+numStrains)=comm_params.MIC_ERT(istrain+totStrains);  
            params_conj_freq(icomm+numStrains)=comm_params.conj_freq(istrain+totStrains);  
            
            if isfield(comm_params, 'colors')
                params_colors(icomm+numStrains,:)=comm_params.colors(istrain,:); 
            else
                params_colors(icomm+numStrains,:)=[.5 .5 .5]; 
            end
            
       end
    end
    params.cs=params_cs;
    params.Vs=params_Vs;
    params.Ks=params_Ks;
    params.species=params_species;
    params.plasmids=params_plasmids;
    params.strains=params_strains;
    params.colors=params_colors;
    params.numStrains=numStrains;
    params.MIC_AMC=params_MIC_AMC;
    params.MIC_MER=params_MIC_MER;
    params.MIC_IMP=params_MIC_IMP;
    params.MIC_ERT=params_MIC_ERT;
    params.conj_freq=params_conj_freq;
    
    