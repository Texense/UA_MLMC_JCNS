% VErrMeanFig.m

% This file is part of UA_MLMC_JCNS, a collection of
% numerical experiments on multilevel Monte Carlo for
% spiking neuron networks.   It accompanies the
% paper "Multilevel Monte Carlo for Cortical Circuit Models,"
% to appear in the Journal of Computational Neuroscience.
% 
% Copyright (C) 2021 by Zhuo-Cheng Xiao
% <zx555@nyu.edu>
% 
% This program is free software; you can redistribute
% it and/or modify it under the terms of the GNU
% General Public License as published by the Free
% Software Foundation; either version 2 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it
% will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General
% Public License along with this program; if not, write
% to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA 02110-1301 USA.

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