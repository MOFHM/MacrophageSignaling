function [significant_peptides_index,significant_peptides_index_upreg,significant_peptides_index_downreg,fc] = plot_generation2_fullMeasure(final_data_p,final_data_p_transformed_full_measured,int_data,header,header_p,p_names_p,position_protein_p,nr_rep,nr_miss_total,Colors,sample_names,nr_sample_groups,dir)

for j = 1:size(final_data_p,2)
    figure
    WindowAPI(gcf,'maximize')
    for i = 1:size(final_data_p{j},2)
        h = histfit(final_data_p{j}(:,i));
        hold on
        h(1).Visible = 'off';
        h(2).Color = Colors(i,:);
    end
    xlabel('Intensity')
    ylabel('Count')
    title(strcat('Distribution of peptides intensities p',num2str(j)))
    grid on
    grid minor
    exportgraphics(gcf,strcat(dir,'/Distribution of peptides intensities_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

for j = 1:size(final_data_p_transformed_full_measured,2)
    figure
    WindowAPI(gcf,'maximize')
    for i = 1:size(final_data_p_transformed_full_measured{j},2)
        h = histfit(final_data_p_transformed_full_measured{j}(:,i));
        hold on
        h(1).Visible = 'off';
        h(2).Color = Colors(i,:);
    end
    xlabel('Intensity')
    ylabel('Count')
    title(strcat('Distribution of peptides intensities p',num2str(j),' zscore'))
    grid on
    grid minor
    exportgraphics(gcf,strcat(dir,'/Distribution of peptides intensities_p',num2str(j),'_zscore','.pdf'),'ContentType','vector')
end
close all

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

for j = 1:size(final_data_p_transformed_full_measured,2)
    corrplot(final_data_p_transformed_full_measured{j},'type','Spearman','tail','both','testR','on','varNames',header_p{j})
    WindowAPI(gcf,'maximize')
    exportgraphics(gcf,strcat(dir,'/Corr_plot_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

for j = 1:size(final_data_p_transformed_full_measured,2)
    ams = clustergram(final_data_p_transformed_full_measured{j},'RowLabels',p_names_p{j},'ColumnLabels',header_p{j},'ImputeFun',@knnimpute,'ColumnLabelsRotate',45,'Standardize','column')
	plot(ams)
    exportgraphics(gcf,strcat(dir,'/Clustergram_plot_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all


fc = cell(size(final_data_p_transformed_full_measured,2),nchoosek(nr_sample_groups,2));
ind_fc = 0;
%qnorm gives better results
for j = 1:size(final_data_p_transformed_full_measured,2)
    ind_fc = 0;
    for i = 1:nr_rep:size(final_data_p_transformed_full_measured{j},2)-nr_rep
        for k = i+nr_rep:nr_rep:size(final_data_p_transformed_full_measured{j},2)
            ind_i = (i-1)/nr_rep+1;
            ind_k = (k-i)/nr_rep+1;
            WindowAPI(gcf,'maximize')
            %With permutations
            % pp = mattest(final_data_p_transformed_full_measured{j}(:,i:i+nr_rep-1),final_data_p_transformed_full_measured{j}(:,k:k+nr_rep-1),'Permute',1000)
            %With unpaired (independent) two tail t-test
            pp = mattest(final_data_p_transformed_full_measured{j}(:,i:i+nr_rep-1),final_data_p_transformed_full_measured{j}(:,k:k+nr_rep-1));
            %Correct for the high false discovery rate - storey method first, BH second
            pp = mafdr(pp); % This method is not centered in 0 for the probability axes
            % pp = mafdr(pp,'BHFDR','true');
            %pp = fdr_BY(pp,0.05,'ind',false)
            mavolcanoplot(final_data_p_transformed_full_measured{j}(:,i:i+nr_rep-1),final_data_p_transformed_full_measured{j}(:,k:k+nr_rep-1),pp,'PlotOnly','on')
            title(strcat('Volcano plot based on t test ',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1}))
            WindowAPI(gcf,'maximize')
            exportgraphics(gcf,strcat(dir,'/Volcano_plot_p',num2str(j),'_zscore_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector')
            ind_fc = ind_fc+1;
            for row = 1:size(final_data_p_transformed_full_measured{j},1)
                vec_comp1 = final_data_p_transformed_full_measured{j}(row,i:i+nr_rep-1);
                vec_comp2 = final_data_p_transformed_full_measured{j}(row,k:k+nr_rep-1);
                vec_comp1(find(isnan(vec_comp1))) = -inf;
                vec_comp2(find(isnan(vec_comp2))) = -inf;
                fc{j,ind_fc} = [fc{j,ind_fc},mean(vec_comp1)-mean(vec_comp2)];
            end
            pp(find(isnan(pp))) = 0;
            figure
            WindowAPI(gcf,'maximize')
            up_reg =  ((pp<0.05) & ((fc{j,ind_fc}>1)'));
            down_reg = ((pp<0.05) & ((fc{j,ind_fc}<-1)'));
            %colsc = ((pp<0.05) & ((fc{j,ind_fc}>1)' | (fc{j,ind_fc}<-1)'));
            colsc = up_reg|down_reg;
            significant_peptides_index{ind_i,ind_k+ind_i-1,j} = find(colsc==1);
            significant_peptides_index_upreg{ind_i,ind_k+ind_i-1,j} = find(up_reg==1);
            significant_peptides_index_downreg{ind_i,ind_k+ind_i-1,j} = find(down_reg==1);
            colsc = repmat(colsc,1,3);
            colsc = double(colsc);
            colsc(find(colsc==0)) = -1;
            colsc(find(up_reg==1),:) = repmat([1 0.55 0],numel(find(up_reg==1)),1);
            colsc(find(down_reg==1),:) = repmat([0.4660 0.6740 0.1880],numel(find(down_reg==1)),1);
            colsc(find(colsc(:,1)==-1),:) = repmat([90 60 100]/255,numel(find(colsc(:,1)==-1)),1);
            % colsc(find(colsc==0)) = -1;
            % colsc(find(colsc==1)) = 10;
            % colsc(find(colsc(:,1)==10),:) = repmat([1 0.55 0],numel(find(colsc(:,1)==10)),1);
            % colsc(find(colsc(:,1)==-1),:) = repmat([90 60 100]/255,numel(find(colsc(:,1)==-1)),1);
            scatter(fc{j,ind_fc},-log10(pp),9,colsc,'*')
            grid on
            grid minor
            yline(-log10(0.05),'-','Threshold p-value(5%)','Color','black','LineWidth',1.5,'LabelVerticalAlignment','middle','FontSize',11)
            xline(-1,'-','FC -1','Color','black','LineWidth',1.5,'LabelHorizontalAlignment','center','FontSize',11)
            xline(1,'-','FC 1','Color','black','LineWidth',1.5,'LabelHorizontalAlignment','center','FontSize',11)
            ylabel('-log_1_0(q-value)')
            xlabel('Log_2(FoldChange)')
            ax = gca;
            ax.XAxis.FontSize = 11;
            ax.XLabel.FontSize = 13;
            %ax.XLabel.FontWeight = 'bold';
            ax.YAxis.FontSize = 11;
            ax.YLabel.FontSize = 13;
            %ax.YLabel.FontWeight = 'bold';
            annotation('textbox','String','Upregulated','Color',[1 0.55 0],'FontSize',13,'LineStyle','none','Position',[0.8 0.7 0.1 0.1])
            annotation('textbox','String','Downregulated','Color',[0.4660 0.6740 0.1880],'FontSize',13,'LineStyle','none','Position',[0.15 0.7 0.1 0.1])
            title(strcat("Volcano plot based on t test ",sample_names{ind_i},"vs ",sample_names{ind_k+ind_i-1}),'FontSize',15)
            exportgraphics(gcf,strcat(dir,'/Self_Volcano_plot_p',num2str(j),'_transfromedData_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector')
            
            fcc1 = fc{j,ind_fc};
            fcc1(find(isinf(fcc1))) = NaN;
            figure
            WindowAPI(gcf,'maximize')
            histfit(fcc1)
            grid on
            grid minor
            xlabel('Log Fold change')
            ylabel('Counts')
            title(strcat('Fold change distribution of data for ',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1}))
            exportgraphics(gcf,strcat(dir,'/Fold_change_distribution_p',num2str(j),'_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector')
            
            
            
            pp = mattest(final_data_p_transformed_full_measured{j}(:,i:i+nr_rep-1),final_data_p_transformed_full_measured{j}(:,k:k+nr_rep-1),'Permute',1000);
            figure
            WindowAPI(gcf,'maximize')
            histogram(pp,20)
            yline(size(final_data_p_transformed_full_measured{j},1)*0.05,'r','p-Value Threshold','LineWidth',1.5)
            title(strcat('Histogram of p-values of non-corrected permutation test (1000) for ',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1}))
            ylabel('Counts')
            xlabel('p-Value')
            grid on
            grid minor
            exportgraphics(gcf,strcat(dir,'/Histogram_p_values_Permutation',num2str(j),'_zscore_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector');
            
            pp = mattest(final_data_p_transformed_full_measured{j}(:,i:i+nr_rep-1),final_data_p_transformed_full_measured{j}(:,k:k+nr_rep-1));
            figure
            WindowAPI(gcf,'maximize')
            histogram(pp,20)
            yline(size(final_data_p_transformed_full_measured{j},1)*0.05,'r','p-Value Threshold','LineWidth',1.5)
            title(strcat('Histogram of p-values of non-corrected two sample t-test for ',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1}))
            ylabel('Counts')
            xlabel('p-Value')
            grid on
            grid minor
            exportgraphics(gcf,strcat(dir,'/Histogram_p_values_t_test',num2str(j),'_zscore_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector')
        end
    end
end
close all

for j = 1:size(final_data_p_transformed_full_measured,2)
    cr = corr(final_data_p_transformed_full_measured{j},'type','Spearman');
    % aa = clustergram(cr,'RowLabels',header_p1,'ColumnLabels',header_p1,'ImputeFun',@knnimpute,'ColumnLabelsRotate',45,'LabelsWithMarkers','true','Cluster','row','Symmetric','false')
    % aa.Colormap = 'hot'
    figure
    WindowAPI(gcf,'maximize')
    img = imagesc(cr);
    colormap(copper)
    %colormap(jet)
    colorbar('EastOutside')
    set(gca,'XTick',1:size(final_data_p_transformed_full_measured{j},2),'XTickLabels',header_p{j})
    xtickangle(45)
    ylabel('Phosphopeptide number')
    xlabel('Experiment ID')
    title('Heatmap for quantifiable phosphopeptides')
    exportgraphics(gcf,strcat(dir,'/Corr_heatmap_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

for j = 1:size(final_data_p,2)
    figure
    WindowAPI(gcf,'maximize')
    % violinplot(int_data)
    %al_goodplot(int_data)
    aa = violinplot(final_data_p{j});
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

for j = 1:size(final_data_p_transformed_full_measured,2)
    figure
    WindowAPI(gcf,'maximize')
    % violinplot(int_data)
    %al_goodplot(int_data)
    aa = violinplot(final_data_p_transformed_full_measured{j});
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
    set(gca,'XTick',1:size(final_data_p_transformed_full_measured{j},2),'XTickLabels',header_p{j})
    xtickangle(45)
    grid on
    grid minor
    exportgraphics(gcf,strcat(dir,'/Violin_transfromedData_p',num2str(j),'.pdf'),'ContentType','vector')
end
close all

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


