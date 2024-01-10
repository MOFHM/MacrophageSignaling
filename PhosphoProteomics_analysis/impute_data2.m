function [data_fin,nr_missing_val] = impute_data2(int_data,nr_rep,nr_sample_groups,penalty_factor,distribution,param,missing_values,position_of_mean_in_param,position_of_std_in_param)
%Function that imputes the missing values of a matrix that has nr_rep
%replicates of nr_sample_groups phenotypes. The values are imputed
%by selecting the values that would make the nr_rep points as close as
%possible to a gaussian curve, keeping a limit with regard to the admisable
%std. The std is estimated from the fully measured data points by the
%function find_dist. When nr_rep-1 points are imputed, a penalty factor
%will be applied.
%int_data - initial data with zeros as missing elements
%nr_rep - number of replicate values
%nr_sample_groups - number of analyzed phenotypes
%penalty_factor - the peanlty applied to nr_rep-1 imputations
%distribution - the desired fitted distribution
%param - the parameters of the desired distribution estimated by the fully
%measured datapoints
%position_of_mean/std_in_param - the index corresponding to the mean/std
%from the param matrix.


%data_fin = ones(size(int_data));
%Number of missing values per phenotype and number of phosphates
nr_missing_val = cell(missing_values,2,size(int_data,2)/nr_rep);
nr_missing_val(:,:,:) = {0};
data_fin = {};
for kk = 1:missing_values
    data_fin{kk} = int_data;
for i = 1:nr_rep:size(int_data,2)
    for j = 1:size(int_data,1)
    vec = int_data(j,i:i+nr_rep-1);
    data_fin{kk}(j,i:i+nr_rep-1) = int_data(j,i:i+nr_rep-1);
    [vec_sort,initial_ind] = sort(vec,'descend');
    ind = find(vec_sort==0);
    if (~isempty(ind))
        
        if (length(ind)==kk)
            %Store the location of the imputed value
            nr_missing_val{kk,1,(i-1)/nr_rep+1} = nr_missing_val{kk,1,(i-1)/nr_rep+1}+1;
            nr_missing_val{kk,2,(i-1)/nr_rep+1} = [nr_missing_val{kk,2,(i-1)/nr_rep+1},j];
            
            for ii = 1:kk
              
            %Mean of the replicates values    
            mean_vector_initial = mean(vec_sort(1:ind(ii)-1));
            %Mean of the means for selected distribution on the specific
            %replicate group
            param_mean = param{1,(i-1)/nr_rep+1}{1,1}(1,position_of_mean_in_param);
            %Mean of the std for selected distribution on the specific
            %replicate group
            param_std =  param{1,(i-1)/nr_rep+1}{1,1}(1,position_of_std_in_param);
            
           % param_mean_range = param{1,(i-1)/nr_rep+1}{1,2}(1,position_of_mean_in_param);
           
            %Std of the std for selected distribution on the specific
            %replicate group
            param_std_range = param{1,(i-1)/nr_rep+1}{1,2}(1,position_of_std_in_param);
            %Number of trials to fit the corresponding distribution using
            %the given datapoints
            nr_trials = 30;
            %Spread the points around the mean value. Here is the difference
            %between inpute_data2 and inpute_data1 as here we first select
            %a point from the left of the curve and then from the right.
            trial_val = linspace(abs(mean_vector_initial+(-1)^ii*param_mean/2),abs(mean_vector_initial+(-1)^(ii+1)*param_std_range),nr_trials);
            limit_values = [param_std-param_std_range;param_std+param_std_range];
            
            aic_value = [];
            parameters_fitted_dist = [];
            
            for trial =  1:nr_trials
                
                %Introduce the trial value in the missing slot
                vec_sort(ind(ii)) = trial_val(trial);
                %Fit the desired distribution and extracts its paramters -
                %following the computation of aic value
                [phat,pci] = mle(vec_sort(1:ind(ii)),'Distribtuion',distribution);
                if numel(phat) == 1
					pdfv= pdf(distribution,vec_sort(1:ind(ii)),phat(1));
				elseif numel(phat) == 2
					pdfv= pdf(distribution,vec_sort(1:ind(ii)),phat(1),phat(2));
				else
					pdfv= pdf(distribution,vec_sort(1:ind(ii)),phat(1),phat(2),phat(3));
			    end
				LL =  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
				aic_value(trial) = 2*numel(phat)- 2*LL;
                parameters_fitted_dist = [parameters_fitted_dist;phat];
                
            end
            %The first criteria for filtering the results is that the std
            %of the fitted distribution to do not be larger or smaller that
            %the limit values of the std extracted from the find_dist
            %algorithm
            ind_std1 = find(parameters_fitted_dist(:,position_of_std_in_param)<limit_values(1));
            ind_std2 = find(parameters_fitted_dist(:,position_of_std_in_param)>limit_values(2));
            
            %Lower AIC is better
            aic_value(ind_std1) = inf;
            aic_value(ind_std2) = inf;
            
            %The second criteria for selecting the desired inputation value
            %is based on the aic criterion - select the one with the lowest
            %aic value
            if (isempty(find(aic_value~=inf))) vec_sort(ind(ii)) = mean_vector_initial+(-1)^ii*param_std_range;
            else
            ind_impute = find(aic_value==min(aic_value));
            end
            vec_sort(ind(ii)) = trial_val(ind_impute(1));
            end
            [~,initial_ind_s] = sort(initial_ind);
            vec = vec_sort(initial_ind_s);
            
            %Penalize the imputation of nr_rep-1 values.
            if (kk==missing_values)
                %Extract all the samples for the respective number of
                %phosphosites
                ind_penalty = fix(i/(nr_rep*nr_sample_groups));
                vec_penalty = int_data(j,ind_penalty*nr_rep*nr_sample_groups+1:(ind_penalty+1)*nr_rep*nr_sample_groups);
                %Cut the values for which the imputation was done
                ind_penalty_cut = mod((i:i+nr_rep-1),(nr_rep*nr_sample_groups));
                %Index of 0 means a real index of nr_rep*nr_sample_groups
                ind_penalty_cut(find(ind_penalty_cut==0)) = nr_rep*nr_sample_groups;
                vec_penalty(ind_penalty_cut) = [];
                %Cut the values of zero
                ind_penalty_cut = find(vec_penalty==0);
                vec_penalty(ind_penalty_cut) = [];
                %Extract the means
                mean_penalty = mean(vec_penalty);
                mean_vec = mean(vec(find(vec~=0)));
                %If the mean of the imputed values is higher that the mean
                %of the rest of the measuremets then use the penalty
                if(mean_vec > mean_penalty)
                    penalty = penalty_factor*(mean_vec-mean_penalty);
                    vec(find(vec~=0)) = vec(find(vec~=0))-penalty;
                end
            end
            
            
            data_fin{kk}(j,i:i+nr_rep-1) = vec;
        end
     end
        
    end

end
end
end