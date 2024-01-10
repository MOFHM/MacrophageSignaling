function [best_tr,test_res,boxcox_param] = find_best_transform(data)

list = {};
list{1} = 'sqrt';
list{2} = 'log2';
list{3} = 'log10';
list{4} = 'asinh';
list{5} = 'pinv';
list{6} = 'quantilenorm';
%list{5} = 'acosh';
%list{6} = 'atanh';

%The statistical test want the data as a row
    data = data(find(data~=0));
    data = data(find(data~=-Inf));
    data = data(find(data~=Inf));
    
    [~,boxcox_param] = boxcox(data);
    %Yeo_johnson transformation needed


data_tr = [];
for i = 1:length(list)
    
data_tr(i,:) = feval(list{i},data);


end
%Inverse done separately
data_tr(i+1,:) = 1./data;
list{i+1} = 'inv';

r = {};
for i = 1:length(list)
    
    data_int = data_tr(i,:);
    data_int = data_int(find(data_int~=0));
    data_int = data_int(find(data_int~=-Inf));
    data_int = data_int(find(data_int~=Inf));
    
r{i} = normalitytest(data_int);

end

test_res = [];

for i = 1:length(r)
    
    test_res = [test_res,r{i}(:,1)];
end

[~,ind] = sort(test_res(1,:));

test_res = test_res(:,ind);
best_tr = list(ind);


end