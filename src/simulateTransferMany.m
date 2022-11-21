function [times, ys, t_end, pf] = simulateTransferMany(params, ic, E)
    
    [numDays, ~]=size(E);

    % Solver parameters
    %options = odeset('NonNegative', 1:length(ic));
    options.RelTol = 1e-12;
    options.AbsTol = 1e-12;   
    
    
    %Numerical simulation
    times=[];
    ys=[];
    t_end=-1;
    
    Ti=params.T;
    for i=1:numDays
    
        this_E=E(i,:);
        
        [this_times, this_y] = ode15s(@(t,x)fMany(t,x, params, this_E),[0,1],ic, options);
        
        times=[times; (i-1)+this_times(2:end)];
        ys=[ys; this_y(2:end,:)];
        
        %if params.numStrains>1
        %    Bpi=sum(this_y(:,2:params.numStrains+1),2);
        %    Bfi=sum(this_y(:,params.numStrains+2:end),2);
        %else
        %    Bpi=this_y(:,2:params.numStrains+1);
        %    Bfi=this_y(:,params.numStrains+2:end);
        %end
      
        
        
        iextinct=find(this_y(end,2:end)<params.extinction_threshold);
        ic=this_y(end,:);
        ic(iextinct+1)=0;
        
        
        
        
        
    end
    
    times=params.T.*times;
    
    B=ys(:,2:end);
    BpT=sum(B(end,1:params.numStrains),2);
    BfT=sum(B(end,params.numStrains+1:end),2);
        
    pf=BpT/(BpT+BfT);  %plasmid fraction
    
end

