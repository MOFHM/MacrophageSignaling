function [list] = remove_cont(data)

ind = [];
if (isempty(data)==0)
for i = 1:length(data)
    aa = contains(data{i},"cont",'IgnoreCase',true);
    if (aa == 1)
        ind = [ind;i];
    end
end
end
list = data;
list(ind) = [];
end
    
