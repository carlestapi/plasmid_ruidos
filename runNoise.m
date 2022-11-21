clc
clear all
close all

run('lib/addpath_recurse');
addpath_recurse('lib/');
addpath_recurse('src/');

addpath('/Users/carlestapi/Library/Application Support/MathWorks/MATLAB Add-Ons/Functions/Entropy');
addpath('/Users/carlestapi/Library/Application Support/MathWorks/MATLAB Add-Ons/Functions/Joint Entropy');
addpath('/Users/carlestapi/Library/Application Support/MathWorks/MATLAB Add-Ons/Functions/Mutual Information');

figPath='figures/';
dataPath='series_temporales/';

%% LOAD MODELING PARAMETERS

B0=1e3;
%ic=[1 B0 0];
ic=[1 B0/2 B0/2];

seg_rate=.1;
Vmax=6e-10;
cost=0.6;

MIC=1.1;

%%

params.conj_rate= 0;
params.S0= 1;
params.T= 24;
params.d= 0.1;
params.cs= [877563800 984410100];
params.Ks= [1 1];
params.species={'E'  'E'};
params.plasmids= {'TC'  'WT'};
params.strains={'X1'  'X1'};
%params.colors= [[0.2510    0.8078    0.8902], [0.2510    0.8078    0.8902]];
params.colors= [[0.25 0.80 0.54], [0.25 0.80 0.54]];
params.numStrains= 1;
params.MIC_AMC= [32768 32];
params.MIC_MER= [4 1];
params.MIC_IMP= [16 0.5000];
params.MIC_ERT= [16 0.5000];
params.conj_freq= [0 0];


dose_max=32768; 
params.MIC_AMC=[dose_max*MIC, 32]; 

params.seg_rate=seg_rate;

params.Vs= [cost*Vmax Vmax];
params.extinction_threshold=1e1;
     
%%

files = dir(fullfile(dataPath, '/**/*.csv'));

%BTs=zeros(CA);
%Pfs=zeros(CA);

% Display the names

%%

I_BpT = zeros(1,length(files));
I_BfT = zeros(1,length(files));
I_Bpf = zeros(1,length(files));
    
H_BpT = zeros(1,length(files));
H_BfT = zeros(1,length(files));

Hurst_BpT=zeros(1,length(files));
Hurst_BpF=zeros(1,length(files));
Hurst_E=zeros(1,length(files));

for inoise=1:length(files)
   this_csv=files(inoise).name;
   path_csv=files(inoise).folder;
   [~, this_noise, ~] = fileparts(this_csv);
   T = readtable([path_csv,'/',this_csv]);
   
   E_sample=(T{1:1:end,1});  %Temp!!!!!!!!!!!!!!!!!!!!!
   
   E=E_sample;
   %min-max renormalization
   %E=(E_sample-min(E_sample))/(max(E_sample)-min(E_sample));
   
   disp([num2str(inoise),') ',this_noise,': ',num2str(E')]);
    
    %Simulate experiment
    [times, ys, t_end, ~] = simulateTransferMany(params, ic, E);
    
    %Estimate plasmid fraction in each transfer
    B=ys(:,2:end);
    
    Ti=find(mod(times,params.T)==0);
   
    daily_BpT=sum(B(Ti,1:params.numStrains),2);
    daily_BfT=sum(B(Ti,params.numStrains+1:end),2);
    
    
    %Pfs(:, inoise)=daily_BpT./(daily_BpT+daily_BfT);  %plasmid fraction
    
    %BTs(:, inoise)=log10(daily_BpT(:)+daily_BfT(:));  %Final bacterial density

    %Plot everything
    %plotManyGrowthCurves(times, ys, params, E); %, max_muKs, max_rhos
    
   

    %pause(1);
%     export_fig([figPath,'series_temporales/',this_noise,'.png']);
%     pause(1);
%   close();
    

%         Auto-Correlation BpT-BfT y E

        [acfBpT,lags] = autocorr(daily_BpT,'NumLags',100,'NumSTD',2);
        [acfBfT,lags] = autocorr(daily_BfT,'NumLags',100,'NumSTD',2);
        [acfE,lags] = autocorr(E,'NumLags',100,'NumSTD',2);

        figure;
        plot(lags,acfBpT,'b',lags,acfBfT,'r',lags,acfE,'g');
        ylim([0.0 1.0])
        xlim([0.0 100])
        legend('ACorr BpT','ACorr BfT','ACorr E');
        shg;
        name = "AutoCorrelation BpT(Blue), BfT(Red), E(Green):" + this_noise;
        title(name);
        xlabel('lag');
        ylabel('AC (Normalized)') ;
        temp=['figures/autocorr/',this_noise,'.png'];
        saveas(gca,temp);
        pause(1);
        close();

    %     %Cross-correlation BpT-BfT y E
%         name = "Cross-Correlation BpT(Blue) / BfT(Red) & E:" + this_noise;
%         [cP,lags] = xcorr(daily_BpT,E,'normalized');
%          %Cross-correlation BfT y E
%         [cF,lags] = xcorr(daily_BfT,E,'normalized');
%          %Cross-correlation BpT y Bft
%         [cPF,lags] = xcorr(daily_BfT,daily_BfT,'normalized');
%     
%     
%         figure;
%         plot(lags,cP,'b',lags,cF,'r',lags,cPF,'g');
%         ylim([0.0 1.0])
%         legend('CC BpT y E','CC BfT y E','CC BpT y Bft');
%         shg;
%         title(name);
%         xlabel('lag');
%         ylabel('CC (Normalized)') ;
%         temp=['figures/crosscorr/',this_noise,'.png'];
%         saveas(gca,temp);
%         pause(1);
%         close();
    
%         
%     %Cross-cov BpT-BfT y E
%     name1 = "Cross-Covariance BpT(Blue) / BfT(Red) & E:" + this_noise;
%     [covP,lags] = xcov(daily_BpT,E,'normalized');
%      %Cross-cov BfT y E
%     [covF,lags] = xcov(daily_BfT,E,'normalized');
%      %Cross-cov BpT y Bft
%     [covPF,lags] = xcov(daily_BfT,daily_BfT,'normalized');
%     
%     figure;
%     plot(lags,covP,'b-',lags,covF,'r-',lags,covPF,'g-');
%     legend('CCOV BpT y E','CCOV BfT y E','CCOV BpT y Bft');
%     shg;
%     title(name1);
%     xlabel('lag');
%     ylabel('CCOV (Normalized)') ;
%     temp=['figures/crosscov/',this_noise,'.png'];
%     saveas(gca,temp);
% %     
%     %Matrix Pearson's correlation coefficients
    
%     X = [daily_BpT, daily_BfT, E];
% %     R = corrplot(X);
%     [R,PValue] = corrplot(X,Type="Pearson",TestR="on")
%     temp=['figures/pcorr/',this_noise,'.png'];
%     saveas(gca,temp);
%     
%     
% %     
    
%     %Entropy 
%     
%     H_BpT(inoise) = Entropy(daily_BfT);
%     H_BfT(inoise) = Entropy(daily_BpT);
%     
%     %Joint Entropy
%    % H_Bpf(inoise) = JointEntropy(daily_BfT,daily_BpT); ?
%     
%     %Mutual Information
%      
%     I_BpT(inoise) = MutualInformation(daily_BpT, E);
%     I_BfT(inoise) = MutualInformation(daily_BfT, E);
%     I_Bpf(inoise) = MutualInformation(daily_BpT, daily_BfT);
% %     
% %     % Generalized Hurst Exponent
% %     
%     Hurst_BpT(inoise) = genhurst(daily_BpT);
%     Hurst_BpF(inoise) = genhurst(daily_BfT);
%     Hurst_E(inoise) = genhurst(E);


    % Power Spectrum
%     
%     
%     sig1 = daily_BpT;
%     sig2 = daily_BfT;
%     sig3 = E;
% 
%     [P1, f1] = pspectrum(sig1);
%     [P2, f2] = pspectrum(sig2);
%     [P3, f3] = pspectrum(sig3);
% 
%     
%     namets = "Time Series: " + this_noise;
%     figure
%     t = (0:numel(sig1)-1);
%     subplot(3,2,1)
%     semilogy(t,sig1,'r');
%     ylabel('BpT')
%     grid on
%     title(namets)
%     subplot(3,2,3)
%     semilogy(t,sig2,'g')
%     ylabel('BfT')
%     grid on
%     subplot(3,2,5)
%     plot(t,sig3,'k');
%     ylabel('E')
%     grid on
%     xlabel('Time (secs)')
%     
%     namepspec = "Power Spectrum: " + this_noise;
%     subplot(3,2,2)
%     loglog(f1,P1,'r')
%     ylabel('P1')
%     grid on
%     axis tight
%     title(namepspec)
%     subplot(3,2,4)
%     loglog(f2,P2,'g')
%     ylabel('P2')
%     grid on
%     axis tight
%     
%     subplot(3,2,6)
%     f3log=log(f3);
%     P3log=log(P3);
%     loglog(f3,P3,'k');
%    
%     
%     ylabel('P3')
%     grid on
%     axis tight
%     xlabel('Frequency (Hz)')
%     temp2=['figures/pspectrum/',this_noise,'.png'];
%     saveas(gca,temp2);
%     pause(1);
%     close();
%     
%     
%     
%     
%     %Coherencia
%     [C1xy,f] = mscohere(sig1,sig2);
%    
%     [C2xy,f] = mscohere(sig1,sig3);
%     
%     [C3xy,f] = mscohere(sig2,sig3);
%     
%     figure
%     namecoh = "Coherence Estimate: " + this_noise;
%     plot(f,C1xy,'k',f,C2xy,'r',f,C3xy,'g')
%     title(namecoh)
%     legend('Coh Bp-Bf','Coh Bp-E','Coh Bf-E')
%     grid on
%     %hgca = gca;
%     temp3=['figures/coherence/',this_noise,'.png'];
%     saveas(gca,temp3);
%     
    
%     pause(1);
%     export_fig([figPath,'series_temporales/',this_noise,'.png']);
%     pause(1);
%     close();
   
    
end



