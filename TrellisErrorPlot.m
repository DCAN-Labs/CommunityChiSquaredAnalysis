function TrellisErrorPlot(observed_marginal_mean,observed_standard_error,output_directory,colors)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[nrows,ncols,ngroups] = size(observed_marginal_mean);
if exist('colors','var') == 0
    colors = repmat([0 0 0],ngroups,1);
end
if exist('output_directory','var') == 0
    output_directory = pwd;
end
if isempty(output_directory)
    output_directory = pwd;
end
h = figure(1);
count = 0;
for i = 1:nrows
    for j = 1:ncols
        count = count + 1;
        if i >= j
            subplot(nrows,ncols,count);
            errorbar(1,observed_marginal_mean(i,j,1),observed_standard_error(i,j,1),'Color',colors(1,:),'LineWidth',6)
            hold;
            for k = 2:ngroups
                errorbar(k,observed_marginal_mean(i,j,k),observed_standard_error(i,j,k),'Color',colors(k,:),'LineWidth',6)
            end
            hold;
            if isnan(observed_marginal_mean(i,j,1)) == 0
                set(gca,'Color',[0.6 0.6 0.6])
            end
        end
    end
end
set(gcf,'Position',[0 0 1920 1080],'PaperUnits','points','PaperPosition',[0 0 1920 1080]);
saveas(h,strcat(output_directory,'/chisquared_trellis.tif'));
end

