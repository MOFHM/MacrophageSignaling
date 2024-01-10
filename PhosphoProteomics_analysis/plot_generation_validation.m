function [] = plot_generation_validation(final_data_p,final_data_p_transform,int_data,header,header_p,p_names_p,position_protein_p,nr_rep,nr_miss_total,Colors,sample_names,nr_sample_groups,dir)
%function that create the quality control plots used to assess the initial
%and normalized/transformed and imputed data.

%Create the heatmaps of the initial data
figure
WindowAPI(gcf,'maximize')
int_data(find(int_data==0)) = nan;
set(groot,'defaultAxesTickLabelInterpreter','none');
img = imagesc(int_data,'AlphaData',~isnan(int_data));
colormap(copper)
%colormap(jet)
colorbar('EastOutside')
set(gca,'XTick',1:size(int_data,2),'XTickLabels',header)
set(gca,'yTick',1:size(int_data,2),'XTickLabels',header)
xtickangle(45)
ylabel('Experiment ID')
xlabel('Experiment ID')
title('Heatmap for quantifiable phosphopeptides')
set(gca,'Color',[170 170 170]/255)
exportgraphics(gcf,strcat(dir,'/Heatmap_initial_data.pdf'),'ContentType','vector')
close all


%Create the correlation plot of the quantile normalized data
for j = 1:size(final_data_p_transform,2)
    corrplot(final_data_p_transform{j},'type','Spearman','tail','both','testR','on','varNames',header_p{j})
    WindowAPI(gcf,'maximize')
    exportgraphics(gcf,strcat(dir,'/Corr_plot_qnorm_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all


%Create the clustergram for the quantile normalized data
for j = 1:size(final_data_p_transform,2)
    aux = final_data_p_transform{j};
    aux(find(isnan(aux))) = 0;
    aux = impute_data_normal(aux);
    ams = clustergram(aux,'RowLabels',p_names_p{j},'ColumnLabels',header_p{j},'ImputeFun',@knnimpute,'ColumnLabelsRotate',45,'Standardize','column')
	plot(ams)
    exportgraphics(gcf,strcat(dir,'/Clustergram_plot_qnorm_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all


%Create the correlation heatmap for the quantile normalized data
for j = 1:size(final_data_p_transform,2)
    aux = final_data_p_transform{j};
    aux = impute_data_normal(aux);
    cr = corr(aux,'type','Spearman');
    % aa = clustergram(cr,'RowLabels',header_p1,'ColumnLabels',header_p1,'ImputeFun',@knnimpute,'ColumnLabelsRotate',45,'LabelsWithMarkers','true','Cluster','row','Symmetric','false')
    % aa.Colormap = 'hot'
    figure
    WindowAPI(gcf,'maximize')
    img = imagesc(cr);
    colormap(copper)
    %colormap(jet)
    colorbar('EastOutside')
    set(gca,'XTick',1:size(final_data_p_transform{j},2),'XTickLabels',header_p{j})
    xtickangle(45)
    ylabel('Phosphopeptide number')
    xlabel('Experiment ID')
    title('Heatmap for quantifiable phosphopeptides')
    exportgraphics(gcf,strcat(dir,'/Corr_heatmap_qnorm_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

%Create the intensity distribution of the initial data
for j = 1:size(final_data_p,2)
    figure
    WindowAPI(gcf,'maximize')
    % violinplot(int_data)
    %al_goodplot(int_data)
    aux = final_data_p{j};
    aux(find(aux==0)) = NaN;
    aa = violinplot(aux);
    for i = 1:length(aa)
        aa(i).ViolinColor = Colors(i,:);
    end
    ylabel('log_2(Intensity)')
    xlabel('Experiment ID')
    title('Intensity distributios of raw data','FontSize',15)
    set(gca,'XTick',1:size(final_data_p{j},2),'XTickLabels',header_p{j})
    xtickangle(45)
    ax = gca;
    ax.XAxis.FontSize = 11;
    ax.XLabel.FontSize = 13;
    ax.XLabel.FontWeight = 'bold';
    ax.YAxis.FontSize = 11;
    ax.YLabel.FontSize = 13;
    ax.YLabel.FontWeight = 'bold';
    grid on
    grid minor
    exportgraphics(gcf,strcat(dir,'/Violin_ini_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

%Create the intensity distribution of the quantile normalized data
for j = 1:size(final_data_p_transform,2)
    figure
    WindowAPI(gcf,'maximize')
    % violinplot(int_data)
    %al_goodplot(int_data)
    aa = violinplot(final_data_p_transform{j});
    for i = 1:length(aa)
        aa(i).ViolinColor = Colors(i,:);
    end
    ylabel('log_2(Intensity)')
    xlabel('Experiment ID')
    title('Intensity distributios of normalized data','FontSize',15)
    ax = gca;
    ax.XAxis.FontSize = 11;
    ax.XLabel.FontSize = 13;
    ax.XLabel.FontWeight = 'bold';
    ax.YAxis.FontSize = 11;
    ax.YLabel.FontSize = 13;
    ax.YLabel.FontWeight = 'bold';
    set(gca,'XTick',1:size(final_data_p_transform{j},2),'XTickLabels',header_p{j})
    xtickangle(45)
    grid on
    grid minor
    exportgraphics(gcf,strcat(dir,'/Violin_qnorm_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

%Extract the number of missing values per phenotype and per nubmer of
%phosphorylated sites
nrm = nr_miss_total(:,2,:);
nrm = squeeze(nrm);
for i = 1:size(nrm,2)
    for j = 1:size(nrm,1)
        nrm{j,i} = numel(nrm{j,i});
    end
end
for i = 1:size(nrm,2)
    figure
    WindowAPI(gcf,'maximize')
    bar([nrm{:,i}])
    grid on
    xlabel('Number of missing values')
    ylabel('Number of phospho peptides')
    set(gca,'XTick',1:5,'XTickLabels',[0,1,2,3,4])
    title(strcat('Missing values exp',num2str(i)))
    exportgraphics(gcf,strcat(dir,'/Number of missing values - exp_',num2str(i),'.pdf'),'ContentType','vector');
end
close all

%Distribution of mean intensities function of missing values
nrm = nr_miss_total(:,2,:);
nrm = squeeze(nrm);
int_data(find(int_data==0)) = nan;
for ii = 1:size(nrm,2)
    aa = {};
    for j = 1:size(nrm,1)
        aa{j} = mean(int_data(nrm{j,ii},ii+(ii-1)*3:ii+ii*3),2,'omitnan');
    end
    figure
    WindowAPI(gcf,'maximize')
    for i = 1:size(aa,2)-1
        h = histfit(aa{:,i});
        hold on
        % h(1).Visible = 'off'
        h(1).FaceAlpha = 0.1;
        h(1).EdgeAlpha= 0.1;
        %h(1).FaceColor = Colors(i,:);
        %h(2).Color = Colors(i,:);
        h(2).Color = h(1).FaceColor;
    end
    xlabel('meanExpression')
    ylabel('count')
    grid on
    grid minor
    l = legend('0','0','1','1','2','2','3','3');
    title(l,'Number of missing values')
    l.Box = 'off';
    title(strcat('Distribution of mean intensities function of missing values for exp',num2str(ii)))
    exportgraphics(gcf,strcat(dir,'/Distribution of missing values_',num2str(ii),'.pdf'),'ContentType','vector')
end
close all

%Average intensity in sample with respect to the total intensity
int_data1 = int_data;
int_data1(find(int_data==0)) = nan;
ss = mean(int_data1,'omitnan')-mean(int_data1(:),'omitnan');
figure
WindowAPI(gcf,'maximize')
bar(ss,'black')
ylabel('sample averae - total average')
xlabel('Experiment ID')
title('Average intensity in sample with respect to the total intensity')
set(gca,'XTick',1:length(header),'XTickLabels',header)
xtickangle(45)
grid on
grid minor
yline(2.5,'-','Thresholod','Color','r','LineWidth',3,'LabelVerticalAlignment','middle')
yline(-2.5,'-','Thresholod','Color','r','LineWidth',3,'LabelVerticalAlignment','middle')
ylim([-3 3])
exportgraphics(gcf,strcat(dir,'/Average_intensity.pdf'),'ContentType','vector')
close all


