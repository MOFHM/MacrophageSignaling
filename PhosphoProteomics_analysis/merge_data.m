function [merged_data] = merge_data(data_fin,missing_values,nr_missing_val,nr_rep)

merged_data = data_fin;

for j = 2:missing_values
    for i = 1:nr_rep:size(data_fin{1,1},2)
        
        k = j;
        while(k>1)    
    for ii = 2:numel(nr_missing_val{k-1,2,(i-1)/nr_rep+1})
        merged_data{1,j}(nr_missing_val{k-1,2,(i-1)/nr_rep+1}(ii),i:i+nr_rep-1) =  merged_data{1,k-1}(nr_missing_val{k-1,2,(i-1)/nr_rep+1}(ii),i:i+nr_rep-1);
    end
    k = k-1;
        end
        
    end
end
end