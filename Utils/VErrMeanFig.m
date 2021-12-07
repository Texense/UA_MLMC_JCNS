%% Ploting Mean+std for Error vs time

function [] = VErrMeanFig(PlotX, YDatas, xlab, ylab, tit, leg, color)
hold on

M1 = mean(YDatas); S1 = std(YDatas);
plot(PlotX,M1,'-','LineWidth',3, 'color',color)
hh1 = plot(PlotX,M1+S1,'--', 'color',color);
hh2 = plot(PlotX,M1-S1,'--', 'color',color);
set(get(get(hh1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(hh2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xlabel(xlab)
ylabel(ylab)
title(tit)
legend(leg,'Interpreter','tex');
hold off
end