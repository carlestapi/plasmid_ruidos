 function plotManyGrowthCurves(time, y, params, E, max_muKs, max_rhos, gm, dyn)


if nargin<4
    strain_labels={};
    strain_colors=jet(params.numStrains);
else
    if isfield(params,'strains')
        strain_labels=params.strains;
    else
        strain_labels={};
    end
    
    if isfield(params,'colors')
        strain_colors=params.colors;
       strain_colors=jet(params.numStrains);
    else
       strain_colors=jet(params.numStrains);
    end
end

if nargin<5
    max_muKs=10e-10;
    max_rhos=12e8;
end


if nargin<7
    gm=[];
end

if nargin<8
   dyn=nan; 
end

%level_set=[1];
level_color=[.8 .8 .8];
rmax=2400;
rmin=24;

[numDays, numDrugs]=size(E);

figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'fontsize',20);
set(gcf, 'Units','normalized','Position',[0. 0. 1 1]);

T=time(end);

S=y(:,1);
B=y(:,2:end);


[~, totStrains]=size(B);
numStrains=totStrains/2;

tot=sum(B,2);  %Total density of bacteria

maxY=max([max(sum(B(:,1:numStrains),2)), max(sum(B(:,numStrains+1:end),2)) ]);
labels={};
%%
% Plot contour of metacommunity
subaxis(7,3,1,1,1,2,'SpacingHoriz',0.075,'PaddingBottom',.05,'MarginTop',0.01);
if ~isempty(gm)
    fcontour(@(x,y)reshape(pdf(gm,[x(:),y(:)]),size(x)),[0 1],'LineWidth',.5,'LineColor',level_color,'MeshDensity',100); hold on;
end

for i=1:numStrains
    
%    if ~strcmp(strain_labels{i},'X')
        
        %if tot(end)>0

            this_color=strain_colors(i,:);

            %Plot parameters w/ circle radius proportional to relative density
            subaxis(7,3,1,1,1,2,'SpacingHoriz',0.075,'PaddingBottom',.05,'MarginTop',0.01);
            
            this_freqTC=B(end,i)./tot(end);
            this_freqWT=B(end,i+params.numStrains)./tot(end);

            %Klebsiella
            %if strcmp(params.species(i),'K')
            %    str_marker='d';
            %else %E coli
                str_marker='o';
            %end

            rmarker_TC=rmin+(rmax-rmin)*this_freqTC;
            rmarker_WT=rmin+(rmax-rmin)*this_freqWT;
            
            %Plasmid-bearing
            %scatter(params.Vs(i)/max_muKs, params.cs(i)/max_rhos,rmarker_WT,this_color,'filled'); hold on;
            s_TC=scatter(params.Vs(i)/max_muKs, params.cs(i)/max_rhos,rmarker_TC,this_color); hold on;
            s_TC.LineWidth = 1;
            s_TC.MarkerEdgeColor = this_color;
            s_TC.MarkerFaceColor = this_color;
            s_TC.MarkerFaceAlpha = 0.5;
            
            %Plasmid-free
            s_WT=scatter(params.Vs(i+numStrains)/max_muKs, params.cs(i+numStrains)/max_rhos, rmarker_WT,this_color,'MarkerFaceAlpha',1); hold on;
            s_WT.LineWidth = 2;
            s_WT.MarkerEdgeColor = this_color;
            s_WT.MarkerFaceColor = [1 1 1];
            s_WT.MarkerFaceAlpha = 0.01;
            

            %if B(1,i)>0  %Plasmid-free at start of experiment
            %    plot(params.Vs(i)/max_muKs, params.cs(i)/max_rhos, '.k','MarkerFaceColor','k','MarkerSize',8); hold on;
            %end
                
            %plot(params.Vs(i+params.numStrains)/max_muKs,params.cs(i+params.numStrains)/max_rhos, str_marker,'LineWidth',2,'MarkerEdgeColor',this_color,'MarkerSize',rmarker_WT); hold on;
            labels{i}=extractBefore(strain_labels{i},'_');

            set(gca,'fontsize',20);
            xlabel('Specific affinity (\mu/K)','FontSize',24);
            ylabel('Cell efficiency (\rho)','FontSize',24);
            nticks=5;
            axis([0 1 0 1]);
            xticks(0:1/nticks:1);
            xticklabels(0:max_muKs/nticks:max_muKs);
            yticks(0:1/nticks:1);
            yticklabels(0:max_rhos/nticks:max_rhos);

            %Plot density CT vs WT
            subaxis(7,3,1,6,3,2,'SpacingHoriz',0.075,'PaddingBottom',0.001,'PaddingTop',.001);
            %semilogy(time, B(:,i),'-', 'Color', strain_colors(i,:)); hold on;
            %semilogy(time, B(:,i+numStrains),'--', 'Color', strain_colors(i,:)); hold on;
            %yticks(10.^linspace(0,9,10));
            
            plot(time, B(:,i),'-', 'Color', strain_colors(i,:)); hold on;
            plot(time, B(:,i+numStrains),'--', 'Color', strain_colors(i,:)); hold on;
            ylim([0 10^9]);
            
            ylabel('Density','fontsize',24);
            axis([0 time(end) 1 2*maxY]);
            xlabel('Time','fontsize',24);
            if T==24
                xticks(0:T/6:T);
            else
                %xticks(0:(params.T*numDays)/10:params.T*numDays);
                %xticklabels(0:(params.T*numDays)/10:numDays*params.T)
            end
            set(gca,'FontSize',20);
       % end
        
 %   end
    
end
%% CONJ_RATE

subaxis(7,3,2,1,1,2,'SpacingHoriz',0.075,'PaddingBottom',.05,'MarginTop',0.01);
if ~isempty(gm)
    fcontour(@(x,y)reshape(pdf(gm,[x(:),y(:)]),size(x)),[0 1],'LineWidth',.5,'LineColor',level_color,'MeshDensity',100); hold on;
end

for i=1:numStrains
    
%    if ~strcmp(strain_labels{i},'X')
        
        if tot(end)>0

            this_color=strain_colors(i,:);

            %Plot parameters w/ circle radius proportional to relative density
            subaxis(7,3,2,1,1,2,'SpacingHoriz',0.075);
            
            
            %Plasmid-free
            s_WT=bar(i, params.conj_freq(i+numStrains), 'FaceColor',this_color); hold on;
            s_WT.LineWidth = 1;
            %s_WT.MarkerFaceColor = [1 1 1];
            %s_WT.MarkerFaceAlpha = 0.01;
            set(gca,'YScale','log')
            set(gca,'fontsize',20);
            ylabel('Conj freq (log)','FontSize',24);
            
            if(~isempty(find(isnan(params.conj_freq)==1, 1)))
                %axis([0 numStrains+1 0 10*max(params.conj_freq+0.1)]);
                yticks([1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1]);
                %xticklabels({'Strains:'})
                %text(i,-(.05*max(params.conj_freq),params.strains{i}),'FontSize',20,'VerticalAlignment','top','HorizontalAlignment','center');
                %axis([0 numStrains+1 0 10*max(params.conj_freq)]);
            end
            
        end
        xticks(1:numStrains)
        xticklabels(strain_labels)
end

%%

subaxis(7,3,1,3,3,2,'SpacingHoriz',0.075,'PaddingBottom',0.001,'PaddingTop',.001);
Bp=sum(B(:,1:numStrains),2);
Bf=sum(B(:,numStrains+1:end),2);
semilogy(time, abs(Bp),'-', 'Color','k'); hold on;
semilogy(time, abs(Bf),'--', 'Color', 'k'); hold on;

ylabel('Density (log)','fontsize',24);
axis([0 time(end) 1 2*maxY]);
%xlabel('Time','fontsize',24);
if T==24
    xticks(0:T/6:T);
else
    xticks(0:params.T:params.T*numDays);
    %xticklabels(0:params.T:numDays*params.T)
    xticklabels({})
end

yticks(10.^linspace(0,9,10));
set(gca,'FontSize',20);
hleg3=legend({'Plasmid-bearing','Plasmid-free'},'Location','SouthEast');
set(hleg3,'FontSize',20);

%%
num_colors=256;
cmap_grey=flipud(gray(num_colors));
subaxis(7,3,1,5,3,1,'SpacingHoriz',0.075,'SpacingVert',.01,'PaddingBottom',0.01,'PaddingTop',.01);
for d=1:numDrugs
   for n=1:numDays
       if E(n,d)>=0
           this_val=E(n,d);
           icolor=round(this_val*(num_colors-1))+1;
           if icolor<1
               icolor=1;
           elseif icolor>255
              icolor=256; 
           end
          this_color=cmap_grey(icolor,:);
       end
       rectangle('Position',[(n-1) (d-1)  1 1],'FaceColor',this_color,'EdgeColor',this_color)
   end
end
yticks(0.5:1:numDrugs-.5)
ylabel('Drug','fontsize',24);
axis([0 numDays 0 numDrugs]);
%xlabel('Time','fontsize',24);
xticks(0:numDays);
xticklabels({})
yticklabels({'AMC','MER','IMP','ERT'});
set(gca,'FontSize',20);

%%
subaxis(7,3,3,1,1,2,'SpacingHoriz',0.075,'PaddingBottom',.075,'MarginTop',0.01);
for d=1:numDrugs
   for n=1:numStrains
       
       %TC
       if d==1
           this_MIC=params.MIC_AMC(n)/32768;
       elseif d==2
           this_MIC=params.MIC_MER(n)/256;
       elseif d==3
           this_MIC=params.MIC_IMP(n)/1024;
       elseif d==4
           this_MIC=params.MIC_ERT(n)/32;
       else
           this_MIC=0;
       end
         
       this_alpha=this_MIC;
       if this_alpha>1
           this_alpha=1;
       end
       this_color=[strain_colors(n,1) strain_colors(n,2) strain_colors(n,3) this_alpha] ;
       
       rectangle('Position',[(n-1) (d-1)  1 1],'FaceColor',this_color)
       
       
       %WT
       if d==1
           this_MIC=params.MIC_AMC(numStrains+n)/32768;
       elseif d==2
           this_MIC=params.MIC_MER(numStrains+n)/256;
       elseif d==3
           this_MIC=params.MIC_IMP(numStrains+n)/1024;
       elseif d==4
           this_MIC=params.MIC_ERT(numStrains+n)/32;
       else
           this_MIC=0;
       end
         
       this_color=[strain_colors(n,1) strain_colors(n,2) strain_colors(n,3) this_MIC] ;
       rectangle('Position',[(numStrains+n-1) (d-1)  1 1],'FaceColor',this_color)
       
   end
end
yticks(0.5:1:numStrains+numDrugs-.5)
ylabel('MIC','fontsize',24); hold on;

plot([numStrains numStrains], [0 numDrugs],'k','LineWidth',3)
xlim([0 2*numStrains])
%xticks(0.5:1:numStrains+0.5)
%xticklabels(strain_labels)
xticks([numStrains/2, 3*numStrains/2])
xticklabels({'TC','WT'})


%axis([0 numDays 0 numDrugs+1]);
%xlabel('Time','fontsize',24);
%xticks(0:numDays);
%xticklabels({})
yticklabels({'AMC','MER','IMP','ERT'});
set(gca,'FontSize',20);


