function [ p_mat,f_mat ] = GenerateAnovaPvals(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nmats = max(size(varargin));
p_mat = zeros(size(varargin{1},1),size(varargin{1},2));
f_mat = p_mat;
for i = 1:nmats
    nsubs(i,1) = size(varargin{i},3);
end
group_index = 1;
groups = zeros(sum(nsubs),1);
for i = 1:nmats
    groups(group_index:sum(nsubs(1:i)),1) = i;
    group_index = group_index + nsubs(i);
end
clear group_index
for i = 1:size(varargin{1},1) - 1
    for j = i + 1:size(varargin{1},2)
        y = zeros(sum(nsubs),1);
        group_index = 1;
        for k = 1:nmats
            y(group_index:sum(nsubs(1:k)),1) = varargin{k}(i,j,:);
            group_index = group_index + nsubs(k);
        end
        [p_mat(i,j),tbl] = anova1(y,groups,'off');
        f_mat(i,j) = tbl{2,5};
        p_mat(j,i) = p_mat(i,j);
        f_mat(j,i) = f_mat(i,j);        
    end
end
        
end