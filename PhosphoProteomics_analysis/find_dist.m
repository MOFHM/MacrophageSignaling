function [dist_names_group,parameter_mean_values_group,dist_performance_ind_group] = find_dist(int_data,nr_rep)
%function that takes an input matrix made of replicate measurements of
%multiple phenotypes and serches for the best distibtuions that describes
%the fully measured data.
%int_data - the initial data matrix with a number unkonown of phenotypes
%nr_rep - the number of repetitions for each phenotype (must be the same
%for all the phenotypes)

rng(123)
%Generate a random fitting in order to extract the distributions names and
%number of parameters
prep = rand(1,4);
aa = fitmethis1(prep,'criterion','AIC');
aa = struct2table(aa);
nr_misses_zeros = zeros(1,size(int_data,2)/nr_rep);
dist_names_group = {};
parameter_mean_values_group = {};
dist_performance_ind_group = {};

for ii = 1:nr_rep:size(int_data,2)
    dist_indices = [];
    dist_parameters = {};
    dist_names = aa.name;
    for j = 1:size(int_data,1)
        %Extract the replicates per conditions
    vec = int_data(j,ii:ii+nr_rep-1);
    %Fit the distributions only on full measured data
    if (isempty(find(vec==0)))
        %Save the number of identified vectors on which a distribution can
        %be fitted
        nr_misses_zeros(1,(ii-1)/nr_rep+1) = nr_misses_zeros(1,(ii-1)/nr_rep+1) + 1;
        %Fit the distributions
        data_fit = fitmethis1(vec,'criterion','AIC');
        data_fit = struct2table(data_fit);
        
        for k = 1:size(dist_names,1)
            %Extract the distributions parameters
            ind = find(strcmp(data_fit.name,dist_names(k,1)));
%             dist_indices_or(k,(ii-1)/nr_rep*size(int_data,1)+j) = ind;
%             dist_parameters_or{k,(ii-1)/nr_rep*size(int_data,1)+j} = data_fit.par(ind);
            dist_indices(k,nr_misses_zeros(1,(ii-1)/nr_rep+1)) = ind;
            dist_parameters{k,nr_misses_zeros(1,(ii-1)/nr_rep+1)} = data_fit.par(ind);
            
        end
        
    end
    end
%Create the performance index based on the total number of hits
dist_performance_ind = {};
for i = 1:size(dist_indices,1)
[~,dist_performance_ind{i}] = find(dist_indices(i,:)==1);
end
%Create the performance index based on the total score
for i = 1:size(dist_indices,1)
dist_performance_ind{2,i} = sum(dist_indices(i,:)');
end

for i = 1:length(dist_performance_ind) a(i) = length(dist_performance_ind{1,i}); end
%Create the performance index based on the total number of hits

ind_best_distribution = find(a==max(a));

parameter_mean_values = {};
%Extract the statistics of the parameters
for i = 1:size(dist_parameters,1)
    nn = [];
    for j = 1:size(dist_parameters,2)
        if (isempty(dist_parameters{i,j})==1) dist_parameters{i,j}{1,1} = Nan; end
        %Uncomment when neccesarry
       % if (isempty(dist_parameters{i,j})==0)
        nn = [nn;dist_parameters{i,j}{1,1}];
        %end
    end
     parameter_mean_values{i,1} = mean(nn);
     parameter_mean_values{i,2} = std(nn);
     parameter_mean_values{i,3} = median(nn);
end

%If neccessary to filter the results
% exc_ind = [];
% for i = 1:size(parameter_mean_values,1)
%     c1 = parameter_mean_values{i,1};
%     c2 = parameter_mean_values{i,2};
%     c3 = c1./c2;
%     if (find((c3<2)>=1)) exc_ind = [exc_ind,i]; end
% end
% 
% dist_names(exc_ind) = [];
% parameter_mean_values(exc_ind,:) = [];
% dist_performance_ind(:,exc_ind) = [];

[~,s_ind] = sort([dist_performance_ind{2,:}]);

dist_names = dist_names(s_ind);
parameter_mean_values = parameter_mean_values(s_ind,:);
dist_performance_ind = dist_performance_ind(:,s_ind);

dist_names_group{(ii-1)/nr_rep+1} = dist_names;
parameter_mean_values_group{(ii-1)/nr_rep+1} = parameter_mean_values;
dist_performance_ind_group{(ii-1)/nr_rep+1} = dist_performance_ind;

end