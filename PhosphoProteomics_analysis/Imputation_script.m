%Extract the indices that corresponds to up/downreulated peptides comparing all the possible combinations from the dataset
%final_data_p - the initial data for this analysis that is already
%transformed, imputed and split by the number of phosphorylated sites
%test_name - the name of the normalization method that should be used
%'qnorm' for quantilenromalization or 'zscore' for zscore transformation
%nr_rep - number of replicates per sample
%sample_names - the names of the samples
%nr_sample_groups - the number of the analyzed phenotypes
%dir_current2 - saving directory
%Colors - used colors for the plots.

%Initializations
plot_1 = {};
plot_2 = {};
fc = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));
p_values_total = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));
baseMean = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));
p_values_noncorrected = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));
stder_total = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));
df_pooled_total = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));
prot_comp_flag = cell(nr_analyzed_phos,nchoosek(nr_sample_groups,2));

pat_ini_gene = 'sp|'+ wildcardPattern + '|';
pat_fin_gene = '_HUMAN' ;

%Analyze the data for the different number of phosphorylated sites
%independently
for j = 1:nr_analyzed_phos
    ind_fc = 0;
    annotation_gene = cellfun(@(x) extractBetween(x,pat_ini_gene,pat_fin_gene),p_names_p{j},'UniformOutput',false);
    annotation_gene = cellfun(@(x) x(1),annotation_gene);
    annotation_residue = cellfun(@(x) strsplit(x,';'), position_protein_p{j}, 'UniformOutput', false);
    annotation_residue = cellfun(@(x) x(1),annotation_residue);
    annotation_vector = cellfun(@(x,y) strcat(x," [",y,"]"),annotation_gene,annotation_residue,'UniformOutput',false);
    
    for i = 1:nr_rep:nr_rep*nr_sample_groups-nr_rep
        for k = i+nr_rep:nr_rep:nr_rep*nr_sample_groups
            close all
            ind_i = (i-1)/nr_rep+1;
            ind_k = (k-i)/nr_rep+1;
            
            
            if (imputation_method == 0)
                data_comp_ttestnorm = [final_data_p{j}(:,i:i+nr_rep-1),final_data_p{j}(:,k:k+nr_rep-1)];
            else
                %Initial data consturction
                data_comp_ttestnorm = int_data(:,(j-1)*nr_rep*nr_sample_groups+1:j*nr_rep*nr_sample_groups);
                ic = [];
                for jk = 1:size(data_comp_ttestnorm,1)
                    vec = data_comp_ttestnorm(jk,:);
                    if (sum(vec==0)==nr_rep*nr_sample_groups) ic = [ic,jk]; end
                end
                data_comp_ttestnorm(ic,:) = [];
                data_comp_ttestnorm = data_comp_ttestnorm(:,[i:i+nr_rep-1 k:k+nr_rep-1]);
            end
            data_comp_ttestnorm(find(data_comp_ttestnorm==0))= NaN;
            WindowAPI(gcf,'maximize')
            
            if(isequal(test_name,'qnorm'))
                data_comp_ttestnorm = quantilenorm(data_comp_ttestnorm);
            end
            if(isequal(test_name,'zscore'))
                data_comp_ttestnorm = mean(std(data_comp_ttestnorm,1,1,'omitnan'),'omitnan').*normalize(data_comp_ttestnorm,1);
            end
            if(isequal(test_name,'median'))
                data_comp_ttestnorm = data_comp_ttestnorm - repmat(median(data_comp_ttestnorm,1,'omitnan'),size(data_comp_ttestnorm,1),1);
            end
            if(isequal(test_name,'mean'))
                data_comp_ttestnorm = data_comp_ttestnorm - repmat(mean(data_comp_ttestnorm,1,'omitnan'),size(data_comp_ttestnorm,1),1);
            end
            
            %Impute the data
            for ik = 1:nr_rep:2*nr_rep
                vec_i = data_comp_ttestnorm(:,ik:(ik+nr_rep-1));
                data_statistics = mean(vec_i,'omitnan');
                data_statistics = [data_statistics; std(vec_i,'omitnan')];
                index_remove = [];
                index_impute = [];
                for imx = 1:size(vec_i,1)
                    vec = vec_i(imx,:);
                    stat = data_statistics;
                    ind = find(isnan(vec));
                    if(isempty(ind)==0)
                        if(length(find(ismember(length(ind),distirbution_imputation_penalty)))>0)
                            for il = 1:length(ind)
                                rng(imx+il+2,'twister')
                                vec(ind(il)) = normrnd(stat(1,ind(il))-Meanshift_imputation_penalty*stat(2,ind(il)),STDshift_imputation_penalty*stat(2,ind(il)));
                            end
                        else if(length(find(ismember(length(ind),distirbution_imputation_normal)))>0)
                                for il = 1:length(ind)
                                   rng(imx+il+3,'twister')
                                    vec(ind(il)) = normrnd(stat(1,ind(il))-Meanshift_imputation_normal*stat(2,ind(il)),STDshift_imputation_normal*stat(2,ind(il)));
                                end
                            else if (length(find(ismember(length(ind),imputation_DreamAI)))>0)
                                    index_impute = [index_impute;imx];
                                else if (length(find(ismember(length(ind),imputation_remove)))>0)
                                        index_remove = [index_remove;imx];
                                    end
                                end
                            end
                            
                        end
                    end
                    vec_i(imx,:) = vec;
                end
                if(length(index_impute)>100)
                    Rclear
                    data_R = vec_i(index_impute,:);
                    delete('data_R.txt')
                    delete('NULL')
                    writematrix(data_R,'data_R.txt')
                    data_R = [];
                    Rinit({'DreamAI','limma','R.matlab'},R_exe_location,R_library_location)
                    Rrun('data_R<-read.table("data_R.txt",sep=",")');
                    Rrun('data_R <- apply(data_R,MARGIN = 2,FUN = as.numeric)');
                    Rrun('set.seed(21)');
                    %Rpush('data_R',data_R)
                    Rrun('imputed_data<-DreamAI(data_R,k=4)');
                    data_R = Rpull('imputed_data');
                    if(isempty(data_R)==0)
                        data_R = data_R.Ensemble;
                    end
                    vec_i(index_impute,:) = data_R;
                    Rclear
                else if(length(index_impute)>0)
                        index_remove = [index_remove,index_impute];
                    end
                end
                data_comp_ttestnorm(:,ik:(ik+nr_rep-1)) = vec_i;
            end
            
            
            %%%%Protein abundance compensation%%%%%
            %%%%Protein_abundance_normalization_moment==1 to scale now.
            if(Protein_abundance_normalization_moment==1)
                [p_names_prot,p_names_LFQmean] = normalize_abundance_afterNormalization(protein_file,nr_rep,nr_sample_groups,p_LFQ_columns,header_name_string);
                for ixk = 1:size(data_comp_ttestnorm,1)
                    names_protein_phospho = split(p_names(ixk),';');
                    for jkx = 1:nr_rep:size(data_comp_ttestnorm,2)
                        div_val = [];
                        ff = fix(((jkx-1)/nr_rep+1)/nr_rep);
                        ind_j = mod((jkx-1)/nr_rep+1,nr_rep)+ff;
                        ind_j(find(ind_j==nr_rep))=1;
                        for kxk = 1:length(names_protein_phospho)
                            aux = find(strcmp(p_names_prot,names_protein_phospho(kxk))==1,1);
                            if(isempty(aux)==0) div_val = [div_val,p_names_LFQmean(aux,ind_j)]; end
                        end
                        if(isempty(div_val)) div_val = median(p_names_LFQmean(:,:)); end
                        div_val_mean = mean(div_val);
                        data_comp_ttestnorm(ixk,jkx:jkx+nr_rep-1) = data_comp_ttestnorm(ixk,jkx:jkx+nr_rep-1)-div_val_mean;
                    end
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            data_comp_ttestnormAUX = data_comp_ttestnorm;
            data_comp_ttestnorm(index_remove,:) = [];
			
            if (statistical_test_used == "t_test")
                %Clasically without considering up/down regluated only
                %Compute the ttest for every individual peptide from the two
                %phenotypes
                pp_ini = mattest(data_comp_ttestnorm(:,1:nr_rep),data_comp_ttestnorm(:,1+nr_rep:end));
                %Correct for the high false discovery rate - storey method first, BH second
                %The storey method uses a bootstrap method to select the lambda
                %value (the value from which the p-value distribution resembles
                %the one of the null hypothesis
                %pp_ini = mafdr(pp_ini); % This method is not centered in 0 for the probability axes
                pp_ini = mafdr(pp_ini,'BHFDR','true');
            end
            
            if (statistical_test_used == "moderate_t_test")
                %R version of the t-test!
                Rclear
                data_R_3 = data_comp_ttestnorm;
                delete('data_R.txt')
                delete('NULL')
                writematrix(data_R_3,'data_R.txt')
                pp_ini = [];
                Rinit({'MKmisc','limma','R.matlab'},R_exe_location,R_library_location)
                Rrun('data_R<-read.table("data_R.txt",sep=",")');
                %Rpush('data_R',data_R)
                Rrun('set.seed(21)');
                Rrun('data_R <- apply(data_R,MARGIN = 2,FUN = as.numeric)')
                Rrun('fdr_values<-mod.t.test(data_R,group=factor(c(rep("group 1", ncol(data_R)/2), rep("group 2", ncol(data_R)/2))))');
                Rrun('pp <- fdr_values$adj.p.value');
                Rrun('pp_noncorrected <- fdr_values$p.value');
                pp_ini = Rpull('pp');
                pp_nc = Rpull('pp_noncorrected');
                Rclear
                delete('data_R.txt')
                delete('NULL')
                %%%%%%
            end
            
            pp = ones(size(data_comp_ttestnormAUX,1),1);
            pp(index_remove) = 0.9988759989;
            pp(setdiff(1:length(pp),find(pp==0.9988759989))') = pp_ini;
            data_comp_ttestnorm1 = ones(size(data_comp_ttestnormAUX));
            data_comp_ttestnorm1(index_remove,:) = data_comp_ttestnormAUX(index_remove,:);
            data_comp_ttestnorm1(setdiff(1:length(pp),find(pp==0.9988759989))',:) = data_comp_ttestnorm;
            data_comp_ttestnorm = data_comp_ttestnorm1;
            
            pp_noncorrected = ones(size(data_comp_ttestnormAUX,1),1);
            pp_noncorrected(index_remove) = 0.9988759989;
            pp_noncorrected(setdiff(1:length(pp_noncorrected),find(pp_noncorrected==0.9988759989))') = pp_nc;
            
            ind_fc = ind_fc+1;
            
            %Compute the fold change considering a fold change of INF when
            %having a NaN comparison
            stder = [];
            df_pooled = [];
            data_fc = data_comp_ttestnorm;
            for row = 1:size(data_fc,1)
                vec_comp1 = data_fc(row,1:nr_rep);
                vec_comp2 = data_fc(row,1+nr_rep:end);
                baseMean{j,ind_fc} = [baseMean{j,ind_fc},mean([mean(nonzeros(vec_comp1),'omitnan'),mean(nonzeros(vec_comp2),'omitnan')],'omitnan')];
                fc{j,ind_fc} = [fc{j,ind_fc},mean(vec_comp1,'omitnan')-mean(vec_comp2,'omitnan')];
                std1 = std(vec_comp1,'omitnan');
                std2 = std(vec_comp2,'omitnan');
                if(std1==0 || std2==0)  stder = [stder,0];  df_pooled = [df_pooled,0];
                else
                    n1 = length(isnan(vec_comp1));
                    n2 = length(isnan(vec_comp2));
                    sp = sqrt((std1^2*(n1-1)+std2^2*(n2-1))/(n1+n2-2));
                    %stder = [stder,sqrt((std1/sqrt(n1-1))^2 + (std2/sqrt(n2-1))^2)];
                    stder = [stder,sp*sqrt(1/n1+1/n2)];
                    %df_pooled = [df_pooled,WelchSatterwhaite(std1,std2,n1,n2)];
                    df_pooled = [df_pooled,n1+n2-2];
                end
            end
            
            stder_total{j,ind_fc} = num2cell(stder);
            stder_total{j,ind_fc} = cellfun(@(x) num2str(round(x,7)),stder_total{j,ind_fc},'UniformOutput',false);
            
            df_pooled_total{j,ind_fc} = num2cell(df_pooled);
            df_pooled_total{j,ind_fc} = cellfun(@(x) num2str(round(x,7)),df_pooled_total{j,ind_fc},'UniformOutput',false);
            
            if (Protein_logfc_compensation_flag==1)
                [fc,pp,~] =  protein_logfc_compensation(protein_logfc_compensation_file_name,sample_names,p_names_p{j},fc,stder,df_pooled,pp,1,j,ind_i,ind_k,ind_fc);
                [fc,pp,~] =  protein_logfc_compensation(protein_logfc_compensation_file_name,sample_names,p_names_p{j},fc,stder,df_pooled,pp,-1,j,ind_i,ind_k,ind_fc);
                prot_comp_flag{j,ind_fc} = -1*zeros(size(pp));
                prot_comp_flag{j,ind_fc} = num2cell(prot_comp_flag{j,ind_fc});
                prot_comp_flag{j,ind_fc} = cellfun(@(x) num2str(round(x,7)),prot_comp_flag{j,ind_fc},'UniformOutput',false);
            end
            
            if (Protein_logfc_compensation_flag==2)
                [~,~,ind1] =  protein_logfc_compensation(protein_logfc_compensation_file_name,sample_names,p_names_p{j},fc,stder,df_pooled,pp,1,j,ind_i,ind_k,ind_fc);
                [~,~,ind2] =  protein_logfc_compensation(protein_logfc_compensation_file_name,sample_names,p_names_p{j},fc,stder,df_pooled,pp,-1,j,ind_i,ind_k,ind_fc);
                prot_comp_flag{j,ind_fc} = zeros(size(pp));
                prot_comp_flag{j,ind_fc}(union(ind1,ind2)) = 1;
                prot_comp_flag{j,ind_fc} = num2cell(prot_comp_flag{j,ind_fc});
                prot_comp_flag{j,ind_fc} = cellfun(@(x) num2str(round(x,7)),prot_comp_flag{j,ind_fc},'UniformOutput',false);
            end
            
            p_values_total{j,ind_fc} = num2cell(pp);
            p_values_total{j,ind_fc} = cellfun(@(x) num2str(round(x,7)),p_values_total{j,ind_fc},'UniformOutput',false);
            
            p_values_noncorrected{j,ind_fc} = num2cell(pp_noncorrected);
            % p_values_noncorrected{j,ind_fc} = cellfun(@(x) num2str(round(x,7)),p_values_noncorrected{j,ind_fc},'UniformOutput',false);
            
            
            %p_values_total{j,ind_fc} = [p_values_total{j,ind_fc}{:}];
            pp(find(isnan(pp))) = 0;
            s = figure();
            s.Position = [179 35 910 610];
            %WindowAPI(gcf,'maximize')
            pause(0.5)
            %Plot the data and extract the indicies of the signinificanlty
            %different peptides
            up_reg =  ((pp<test_p_value) & ((fc{j,ind_fc}>test_fc_value)'));
            down_reg = ((pp<test_p_value) & ((fc{j,ind_fc}<-test_fc_value)'));
            %colsc = ((pp<test_p_value) & ((fc{j,ind_fc}>test_fc_value)' | (fc{j,ind_fc}<-test_fc_value)'));
            colsc = up_reg|down_reg;
            significant_peptides_index{ind_i,ind_k+ind_i-1,j} = find(colsc==1);
            significant_peptides_index_upreg{ind_i,ind_k+ind_i-1,j} = find(up_reg==1);
            significant_peptides_index_downreg{ind_i,ind_k+ind_i-1,j} = find(down_reg==1);
            colsc = repmat(colsc,1,3);
            colsc = double(colsc);
            colsc(find(colsc==0)) = -1;
            %colsc(find(up_reg==1),:) = repmat([0.4660 0.6740 0.1880],numel(find(up_reg==1)),1);
            colsc(find(up_reg==1),:) = repmat([0 0 138]/255,numel(find(up_reg==1)),1);
            %colsc(find(down_reg==1),:) = repmat([1 0.55 0],numel(find(down_reg==1)),1);
            colsc(find(down_reg==1),:) = repmat([0 100 0]/255,numel(find(down_reg==1)),1);
            colsc(find(colsc(:,1)==-1),:) = repmat([90 60 100]/255,numel(find(colsc(:,1)==-1)),1);
            scatter(-fc{j,ind_fc},-log10(pp),17,colsc,'*')
            grid on
            grid minor
            yline(-log10(0.05),'-','{\bfFDR} (5%)','Color','black','LineWidth',1.5,'LabelVerticalAlignment','middle','FontSize',13)
            xline(-1,'-','\bfFC -1','Color','black','LineWidth',2.5,'LabelHorizontalAlignment','center','FontSize',15)
            xline(1,'-','\bfFC 1','Color','black','LineWidth',2.5,'LabelHorizontalAlignment','center','FontSize',15)
            %ylabel('-Log_1_0({\itq}-value)')
            ylabel('-Log_1_0(FDR)')
            xlabel('Log_2(FoldChange)')
            ax = gca;
            xlim(max(abs(ax.XLim)).*[-1 1])
            ax.XAxis.FontSize = 17;
            ax.XLabel.FontSize = 21;
            %ax.XLabel.FontWeight = 'bold';
            ax.YAxis.FontSize = 17;
            ax.YLabel.FontSize = 21;
            %ax.YLabel.FontWeight = 'bold';
            %annotation('textbox','String','Upregulated','Color',[0.4660 0.6740 0.1880],'FontSize',17,'LineStyle','none','Position',[0.75 0.7 0.1 0.1])
            %annotation('textbox','String','Downregulated','Color',[1 0.55 0],'FontSize',17,'LineStyle','none','Position',[0.15 0.7 0.1 0.1])
            annotation('textbox','String','Upregulated','Color',[0 100 0]/255,'FontSize',17,'LineStyle','none','Position',[0.75 0.7 0.1 0.1])
            annotation('textbox','String','Downregulated','Color',[0 0 138]/255,'FontSize',17,'LineStyle','none','Position',[0.15 0.7 0.1 0.1])
			title(strcat("Volcano plot based on moderate t test ",sample_names{ind_k+ind_i-1},"vs ",sample_names{ind_i}),'FontSize',17)
			exportgraphics(gcf,strcat(dir_current2,'/Self_Volcano_plot_p',num2str(j),'_',test_name,'_',sample_names{ind_k+ind_i-1},'vs',sample_names{ind_i},'.pdf'),'ContentType','vector')
			
			
			
			

           % vec_search = ["PML","RSAD2","RIPK2","PD1L1","JAK2","JAK3","RELB","NFKB2","SIN3A","JUN","PDPK1","MP2K2","CCR1","NCOA3","AKTS1","AKT2","DPTOR","GAS7","NCOR2","PAK2"];
            %vec_search = ["VIME","STAT1","RSAD2","PML","RIPK2","PD1L1","JAK2","SP100","CD97","XIRP1","TRAF1","NFKB2","SIN3A","ACINU","SRRM2","JUN","JUNB","CCR1","AKT2","GAS7","PAK2","AKT2","BABA1","C163B","SO2B1","NISCH","ITA4"];
            %ind_ann = cellfun(@(x) find(contains(p_names_p{j},x)),vec_search,'UniformOutput',false);
           % ind_ann = cellfun(@(x) x(find(abs(fc{j,ind_fc}(x))==max(abs(fc{j,ind_fc}(x))))),ind_ann);
             add_annotations_to_figure(s,-fc{j,ind_fc},-log10(pp),annotation_vector,20)
            
                        %add_annotations_to_figure(s,-fc{j,ind_fc}(ind_ann),-log10(pp(ind_ann)),annotation_vector(ind_ann),length(ind_ann))

            saveas(gcf,strcat(dir_current2,'/Self_Volcano_plot_p',num2str(j),'_',test_name,'_',sample_names{ind_k+ind_i-1},'vs',sample_names{ind_i},'.fig'))
            
                         %Plot the fold change distribution
                         fcc1 = fc{j,ind_fc};
                         fcc1(find(isinf(fcc1))) = NaN;
                         figure
                         WindowAPI(gcf,'maximize')
                         histfit(fcc1)
                         grid on
                         grid minor
                         xlabel('Log Fold change')
                         ylabel('Counts')
                         ax = gca;
                         ax.XAxis.FontSize = 17;
                         ax.XLabel.FontSize = 21;
                         %ax.XLabel.FontWeight = 'bold';
                         ax.YAxis.FontSize = 17;
                         ax.YLabel.FontSize = 21;
                         %ax.YLabel.FontWeight = 'bold';
                         title(strcat('Fold change distribution of data for ',sample_names{ind_i},"vs ",sample_names{ind_k+ind_i-1}),'FontSize',21)
                         exportgraphics(gcf,strcat(dir_current2,'/Fold_change_distribution_p',num2str(j),'_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector')
            
            
                         %Plot the distribution of p-values if a permutation test is
                         %used
                         pp = mattest(data_comp_ttestnorm(:,1:nr_rep),data_comp_ttestnorm(:,1+nr_rep:end),'Permute',1000);
                         figure
                         WindowAPI(gcf,'maximize')
                         histogram(pp,20)
                         yline(size(data_comp_ttestnorm,1)*0.05,'r','p-Value Threshold','LineWidth',1.5,'FontSize',17)
                         title(strcat("Histogram of p-values of non-corrected permutation test (1000) for ",sample_names{ind_i},"vs ",sample_names{ind_k+ind_i-1}),'FontSize',21)
                         ylabel('Counts')
                         xlabel('p-Value')
                         grid on
                         grid minor
                         ax = gca;
                        ax.XAxis.FontSize = 17;
                         ax.XLabel.FontSize = 21;
                         %ax.XLabel.FontWeight = 'bold';
                         ax.YAxis.FontSize = 17;
                         ax.YLabel.FontSize = 21;
                         %ax.YLabel.FontWeight = 'bold';
                         exportgraphics(gcf,strcat(dir_current2,'/Histogram_p_values_Permutation',num2str(j),'_',test_name,'_',sample_names{ind_i},"vs ",sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector');

                         %Plot the distribution of p-values for the normal t-test
                        pp = mattest(data_comp_ttestnorm(:,1:nr_rep),data_comp_ttestnorm(:,1+nr_rep:end));
                        figure
                         WindowAPI(gcf,'maximize')
                         histogram(pp,20)
                         yline(size(data_comp_ttestnorm,1)*0.05,'r','p-Value Threshold','LineWidth',1.5,'FontSize',17)
                        title(strcat("Histogram of p-values of non-corrected two sample t-test for ",sample_names{ind_i},"vs ",sample_names{ind_k+ind_i-1}),'FontSize',21)
                        ylabel('Counts')
                         xlabel('p-Value')
                         grid on
                         grid minor
                         ax = gca;
                         ax.XAxis.FontSize = 17;
                         ax.XLabel.FontSize = 21;
                         %ax.XLabel.FontWeight = 'bold';
                         ax.YAxis.FontSize = 17;
                         ax.YLabel.FontSize = 21;
                         %ax.YLabel.FontWeight = 'bold';
                         exportgraphics(gcf,strcat(dir_current2,'/Histogram_p_values_t_test',num2str(j),'_',test_name,'_',sample_names{ind_i},'vs',sample_names{ind_k+ind_i-1},'.pdf'),'ContentType','vector')
        end
    end
end
close all

%Plot the distribution of peptide intensities considering the two analyzed
%normalization methods

for j = 1:size(plot_1,2)
    figure
    WindowAPI(gcf,'maximize')
    yt = {};
    xt = {};
    for i = 1:size(final_data_p_zscore{j},2)
        a = final_data_p_zscore{j}(:,i);
        a = a(find(~isnan(a)));
        h = histfit(a);
        hold on
        h(1).Visible = 'off';
        h(2).Color = Colors(i,:);
        %If the probability density is desired
        yt{i} = h(2).YData./sum(h(1).YData);
        xt{i} = h(2).XData;
        %plot(h(2).XData, h(2).YData/sum(yt), 'r', 'LineWidth',2)
        
    end
    xlabel('Intensity')
    ylabel('Count')
    title(strcat('Distribution of peptides intensities p',num2str(j),' zscore'))
    grid on
    grid minor
    exportgraphics(gcf,strcat(dir,'/Distribution of peptides intensities_p',num2str(j),'_zscore','.pdf'),'ContentType','vector')
    figure
    WindowAPI(gcf,'maximize')
    %If the probability density is desired
%     for i = 1:size(final_data_p_zscore{j},2)
%         plot(xt{i},yt{i},'Color',Colors(i,:),'LineWidth',2)
%         hold on
%     end
        
end
close all