function [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal)

statistics = {};
for j = 1:size(score_kin_subExternal,1)
    if(isempty(find(ismember(score_kin_subBackground(:,1),score_kin_subExternal(j,1))))==0)
        chi_kinase = score_kin_subExternal{find(ismember(score_kin_subExternal(:,1),score_kin_subExternal(j,1))),3};
        chi_background = score_kin_subBackground{find(ismember(score_kin_subBackground(:,1),score_kin_subExternal(j,1))),3};
        chi_kinase_total = round(score_kin_subExternal{1,3}/score_kin_subExternal{1,4});
        chi_background_total = round(score_kin_subBackground{1,3}/score_kin_subBackground{1,4});
        
        %Pooled estimate of proportions:
        p0 = (chi_kinase+chi_background) / (chi_kinase_total+chi_background_total);
        
        %Expected counts under the hypothesis of H0:
        chi_kinase_H0 = chi_kinase_total * p0;
        chi_background_H0 = chi_background_total * p0;
        
        observed = [chi_kinase, chi_kinase_total-chi_kinase, chi_background, chi_background_total-chi_background];
        expected = [chi_kinase_H0, chi_kinase_total-chi_kinase_H0, chi_background_H0, chi_background_total-chi_background_H0];
        
        chiScore = sum((observed-expected).^2 ./ expected);
        p_value_chi = 1-chi2cdf(chiScore,1);
        statistics{j,1} = chiScore;
        statistics{j,2} = p_value_chi;
        
        
        as1 = [repmat('a',chi_kinase_total,1); repmat('b',chi_background_total,1)];
        as2 = [repmat(1,chi_kinase,1); repmat(2,chi_kinase_total-chi_kinase,1); repmat(1,chi_background,1); repmat(2,chi_background_total-chi_background,1)];
        
        [tbl,chi2,pval_tab,labels] = crosstab(as1,as2);
        [h,p_fisher,stats] = fishertest(tbl);
        
        statistics{j,3} = stats.OddsRatio;
        statistics{j,4} = p_fisher;
        
    end
    
end

if(isempty(statistics)==0)
statistics_aux = {};
statistics_aux(:,1) = statistics(:,1);
statistics_aux(:,2) = statistics(:,2);
ind = find(cellfun(@isempty,statistics(:,2)));
%Use Storey or BH multiple correction
s1 = num2cell(mafdr([statistics{:,2}]','BHFDR','true'));
%try
%    s1 = num2cell(mafdr([statistics{:,2}]'));
%catch
%    s1 = num2cell(mafdr([statistics{:,2}]','BHFDR','true'));
%end
for k1 = 1:length(ind)
    s1 = [s1(1:ind(k1)-1);'NaN';s1(ind(k1):end)];
end
statistics_aux = [statistics_aux,s1];
statistics_aux(:,4) = statistics(:,3);
statistics_aux(:,5) = statistics(:,4);
ind = find(cellfun(@isempty,statistics(:,4)));
%Use Storey or BH multiple correction
s1 = num2cell(mafdr([statistics{:,4}]','BHFDR','true'));
%try
%    s1 = num2cell(mafdr([statistics{:,4}]'));
%catch
%    s1 = num2cell(mafdr([statistics{:,4}]','BHFDR','true'));
%end
for k1 = 1:length(ind)
    s1 = [s1(1:ind(k1)-1);'NaN';s1(ind(k1):end)];
end
statistics_aux = [statistics_aux,s1];

statistics = statistics_aux;

end