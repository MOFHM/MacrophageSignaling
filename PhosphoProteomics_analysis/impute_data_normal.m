function [data_out] = impute_data_normal(data_in)
%Function that imputes the missing values of a matrix column wise. The
%imputed values are draw from a gaussian distribution shifted in mean by 1.8x the
%std of the initial data distribution (one column at a time) with a new std
%of 0.3x times the initial std.
%The missing values should be either NaNs or zeros!
rng(123)

for i = 1:size(data_in,2)
    %Remove the missing values
    ind_nan = find(isnan(data_in));
    if(isempty(ind_nan))
        vec = data_in(find(data_in(:,i)~=0),i);
    else
        vec = data_in(find(~isnan(data_in(:,i))),i);
    end
    %Fit the distributions
    distribution_columns = fitmethis1(vec,'criterion','AIC');
    %Extract the paramters of the gaussian distribution
    ind = find(strcmp({distribution_columns.name}, 'normal'));
    par = {distribution_columns.par};
    mean_c_d = par{ind}(1);
    std_c_d = par{ind}(2);
    %Impute the new values from the newly created distribution
    for j = 1:size(data_in,1)
        if(data_in(j,i)==0) data_in(j,i) = normrnd(mean_c_d-1.8*std_c_d,0.3*std_c_d); end
    end
    for j = 1:size(data_in,1)
        if(isnan(data_in(j,i))) data_in(j,i) = normrnd(mean_c_d-1.8*std_c_d,0.3*std_c_d); end
    end
    
end

data_out = data_in;



end