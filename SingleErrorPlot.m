function SingleErrorPlot(observed_marginal_mean,observed_standard_error,filename,colors,group_names)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

ngroups = length(observed_marginal_mean);
if exist('colors','var') == 0
    colors = repmat([0 0 0],ngroups,1);
end
h = figure(1);
errorbar(10/100,observed_marginal_mean(1),observed_standard_error(1)*2,'Color',colors(1,:),'LineWidth',10)
fig = gcf;
hold;
fig.InvertHardcopy = 'off';
scatter(10/100,observed_marginal_mean(1),1000,colors(1,:),'filled','MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:));
for k = 2:ngroups
    errorbar((10+(k-1)/ngroups)/100,observed_marginal_mean(k),observed_standard_error(k)*2,'Color',colors(k,:),'LineWidth',10)
    scatter((10+(k-1)/ngroups)/100,observed_marginal_mean(k),1000,colors(k,:),'filled','MarkerEdgeColor',colors(k,:),'MarkerFaceColor',colors(k,:));
end
xlim([0.08 0.13])
set(gca,'Color',[1 1 1]);
hold;
set(gca,'FontName','Arial','FontSize',30)
set(gca,'XTick',(10:1/ngroups:10.9)./100)
if exist('group_names','var') == 0
    group_names = cell(ngroups,1);
    for i = 1:ngroups
        group_names{i} = strcat('group #',num2str(i));
    end
end
set(gca,'XTickLabel',group_names);
set(gca,'XTickLabelRotation',90);
set(gcf,'Position',[0 0 640 480],'PaperUnits','points','PaperPosition',[0 0 640 480]);
saveas(h,strcat(filename,'.tif'));
end

