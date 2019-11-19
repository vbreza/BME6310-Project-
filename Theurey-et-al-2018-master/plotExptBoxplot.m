function plotExptBoxplot
% To plot boxplots of experimental data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Plot Fig 1C for Theurey et al Aging Cell
% Paste relevant data from "Pre-generated simulation data.xls" into x variable
% x = data (col 1), group (col 2) and scatter (col 3)
% Simulations generated using Beard_NC_simulatePopulation.m

% Specify data variable
%data = 'TMRM_Oligo';
%data = 'TMRM_AA';
%data = 'TMRM_Rot';
x = xlsread('Pre-generated data (Fig 1C,2F).xlsx');
data = 'NADH_Rot';   

% Position = figure dimensions for publication
figure('Position', [100 100 120 130])
% Change default font type to Arial
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

% Assign data to y variable (so you don't have to edit code)
y = x;
colours = ('rk');   % Plot expt data in red (changed in Figure to grey), sims in black 

boxplot(y(:,1),y(:,2),'Color',colours,'Symbol','+')
hold on
markersize = 2;
% Plot scatter points next to boxplot
for i = 1:size(colours,2)
    scatter(y(y(:,2)==i,3),y(y(:,2)==i,1),markersize,'filled',colours(i))
end

% Add ylabel and y axis limits
if strcmp(data, 'TMRM_Rot')
    ylabel('\Delta\Psi_m response to Rot. (mV)')
    ylim([90 145])
elseif strcmp(data, 'TMRM_AA')
    ylabel('\Delta\Psi_m response to AA (mV)')
    ylim([90 145])
elseif strcmp(data, 'NADH_Rot')
    ylabel('NAD(P)H response to Rot. (FC)')
    %ylim([90 145])
elseif strcmp(data, 'TMRM_Oligo')
    ylabel('\Delta\Psi_m response to Oligo. (mV)')
    ylim([90 170])
end

% Set x axis limits
xlim([0.5 2.5])

% Set x-axis tick labels
xticklabels({'CNs','Sims'})

% Remove box around plot
box off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%
% % Plot Fig 2F for Aging Cell Revisions
% % x = data (col 1), group (col 2) and scatter (col 3)
% % data from Fig 2F tab in 'Compare 'tg' sims vs expts_AgingCellRev.xls'
% 
% % Position = figure dimensions for publication
% figure('Position', [100 100 220 190])
% % Change default font type to Arial
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultTextFontName', 'Arial')
% 
% % Assign your data to y variable (so you don't have to edit code)
% y = x;
% colours = ('krkrkrkr');   % Plot WTsim in black, DH defect in red 
% 
% boxplot(y(:,1),y(:,2),'Color',colours,'Symbol','+')
% hold on
% markersize = 2;
% % Plot scatter points next to boxplot
% for i = 1:size(colours,2)
%     scatter(y(y(:,2)==i,3),y(y(:,2)==i,1),markersize,'filled',colours(i))
% end
% 
% % Plot means
% scatter(1,nanmean(x(1:100,1)),2,'k')
% scatter(2,nanmean(x(101:200,1)),2,'r')
% scatter(3,nanmean(x(201:300,1)),2,'k')
% scatter(4,nanmean(x(301:400,1)),2,'r')
% scatter(5,nanmean(x(401:500,1)),2,'k')
% scatter(6,nanmean(x(501:600,1)),2,'r')
% scatter(7,nanmean(x(601:700,1)),2,'k')
% scatter(8,nanmean(x(701:800,1)),2,'r')
% 
% % Add ylabel and y axis limits
%     ylabel('Simulated OCR (pmol/O_2/min/\mug protein)')
%     %ylim([90 170])
% 
% % Set x axis limits
% %xlim([0.5 2.5])
% 
% % Set x-axis tick labels
% xticklabels({'Basal','','H+ leak','','ATP synth','','Maximal',''})
% 
% % Remove box around plot
% box off

%%

end
