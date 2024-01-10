function [CV] = CV_comp(final_data_p,nr_rep,Colors,name,method,sample_names,dir_current)

data_cv = final_data_p;
CV = cell(1,size(final_data_p,2));
for i = 1:size(final_data_p,2)
    for j = 1:nr_rep:size(final_data_p{i},2)
        if(method == 1)
cc1 = std(final_data_p{i}(:,j:j+nr_rep-1),0,2,'omitnan')
        else
            cc1 = std(final_data_p{i}(:,j:j+nr_rep-1),0,2,'omitnan')./mean(final_data_p{i}(:,j:j+nr_rep-1),2,'omitnan');
        end
CV{i} = [CV{i},cc1];
    end
    if (method==1)
    cc1 = std(final_data_p{i}(:,1:end),0,2,'omitnan')
    else
        cc1 = std(final_data_p{i}(:,1:end),0,2,'omitnan')./mean(final_data_p{i}(:,1:end),2,'omitnan');
    end
    CV{i} = [CV{i},cc1];
end

for j = 1:size(CV,2)
figure
WindowAPI(gcf,'maximize')
aa = violinplot(CV{j})
for i = 1:length(aa)
aa(i).ViolinColor = Colors(i,:);
end
ylabel(strcat('CV (',name,')'))
xlabel('Expressed phenotype')
title(['Violin plots of CV for ',name])
set(gca,'XTick',1:size(CV{j},2),'XTickLabels',{sample_names{:},'All'})
set(gca,'TickLabelInterpreter','tex')
%xtickangle(45)
grid on
grid minor
saveas(gcf,strcat(dir_current,'/','Violin_CV_',name,num2str(j),'.tif'))
end
close all