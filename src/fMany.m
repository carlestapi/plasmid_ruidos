function fout = fMany(~, y, params, E)

    
    MIC_max=[32768, 256, 1024, 32]; %tmp

    S=y(1);
    numStrains=(length(y)-1)/2;

    B_TC=y(2:numStrains+1);
    B_WT=y(numStrains+2:end);

    if nargin<4
        E=[0];
    end
    
    numDrugs=length(E);
    %disp(['E=',num2str(E)]);
    
    uStot=0;
    dB=zeros(1, numStrains);

    %For plasmid-bearing (TC)
    for i=1:numStrains
        
        uSi=uS(S, params.Vs(i), params.Ks(i));
        uStot=uStot+uSi*B_TC(i);  %Resource consumed

        conj_rate=params.conj_freq(i+numStrains);
        if isnan(conj_rate) %??????????????????
            conj_rate=1e-12; %0; 
        end
        
        %Plasmid transfer B_TC_j -> B_WT_i
        dBconj=0;
        for j=1:numStrains
            dBconj=dBconj+conj_rate*B_TC(j)*B_WT(i);
        end
        
        kill_rate_WT=1-[params.MIC_AMC(i), params.MIC_MER(i), params.MIC_IMP(i), params.MIC_ERT(i)]./MIC_max;
        dBkill=0;
        for d=1:numDrugs
            dBkill=dBkill+(kill_rate_WT(d)*E(d)*B_TC(i));
        end
        dB(i)=params.cs(i)*uSi*B_TC(i)*(1 - params.seg_rate) + dBconj - dBkill - params.d*B_TC(i);
   

        %For plasmid-free (WT)
        uSi=uS(S, params.Vs(numStrains+i), params.Ks(numStrains+i));
        uStot=uStot+uSi*B_WT(i);  %Resource consumed

        conj_rate=params.conj_freq(i+numStrains);
        if isnan(conj_rate) %??????????????????
            conj_rate=0; 
        end
        
        %Plasmid transfer B_TC_j -> B_WT_i
        dBconj=0;
        for j=1:numStrains
                dBconj=dBconj+conj_rate*B_TC(j)*B_WT(i);
        end
        
        kill_rate_TC=1-[params.MIC_AMC(numStrains+i), params.MIC_MER(numStrains+i), params.MIC_IMP(numStrains+i), params.MIC_ERT(numStrains+i)]./MIC_max;
        dBkill=0;
        for d=1:numDrugs
            dBkill=dBkill+(kill_rate_TC(d)*E(d)*B_WT(i));
        end
        
        dB(numStrains+i)=params.cs(numStrains+i)*uSi*B_WT(i) - dBconj - dBkill + params.cs(i)*uSi*B_TC(i)*(params.seg_rate) - params.d*B_WT(i);
        
    end

    dS=-uStot + params.d*(params.S0-S);
    fout = params.T*([dS dB]');
end


function ret = uS(S,V,K)
    ret=(S*V)/(K+S);
end
