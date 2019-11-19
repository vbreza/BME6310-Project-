function plotTimeSeries(tt, stepsize, numSims,...
    plotStateVar, plotStateVar_Choice,stateVarAll,...
    plotStateVarFC, plotStateVarFC_Choice, stateVarFCAll,...
    plotOutput, plotOutput_Choice, OutputAll,...
    plotOutputFC, plotOutputFC_Choice, OutputFCAll)
%%% Plots time series of state variable and outputs
%%% Only plots SVs and Os as chosen by user in calling function

%%% Calculate mean of outputs
if size(stateVarAll,3) == 1 % if only 1 simulation, don't calculate mean!
    meanStateVar = 0;
    meanStateVarFC = 0;
    meanOutput =0;
    meanOutputFC = 0;
else
    meanStateVar = nanmean(stateVarAll,3);
    meanStateVarFC = nanmean(stateVarFCAll,3);
    meanOutput = nanmean(OutputAll,3);
    meanOutputFC = nanmean(OutputFCAll,3);
end

% Get legends for graph titles
[statevar_legend, output_legend] = getLegends;

%plot_time = inf;        % To constrain graph x-axes; inf sets it to auto
plot_cols = 4;

if plotStateVar == 1
    plot_rows = ceil(size(plotStateVar_Choice,2)/plot_cols);
    figure('Name','State vars timeseries (abs)')
    for j = 1:size(plotStateVar_Choice,2)
      plotIndex = plotStateVar_Choice(j);
      subplot(plot_rows,plot_cols,j), hold on
      for i = 1:numSims
        plot(tt,stateVarAll(:,plotIndex,i)); 
      end
      if meanStateVar == 0
      else
        plot(tt,meanStateVar(:,plotIndex),'k','LineWidth',2);
      end
      title(statevar_legend(plotIndex))
      %axis([-inf plot_time -inf inf])
    end
    %Plot deltaPsim 
    figure, hold on
    for i = 1:numSims, plot(tt,stateVarAll(:,19,i)); end
    if meanStateVar == 0
    else
      plot(tt,meanStateVar(:,19),'k','LineWidth',2);
    end
    axis([0 65 0 180])
    xlabel('Time (min)')
    ylabel('Mito. memb. potential (mV)')
end

if plotStateVarFC == 1
    plot_rows = ceil(size(plotStateVarFC_Choice,2)/plot_cols);
    figure('Name','State vars timeseries (foldchange)')
    for j = 1:size(plotStateVarFC_Choice,2)
      plotIndex = plotStateVarFC_Choice(j);
      subplot(plot_rows,plot_cols,j), hold on
      plot(tt,squeeze(stateVarFCAll(:,plotIndex,:))); 
      if meanStateVarFC == 0
      else
        plot(tt,meanStateVarFC(:,plotIndex),'k','LineWidth',2);
      end
      title(statevar_legend(plotIndex))
      %axis([-inf plot_time -inf inf])
    end
    %Plot deltaPsim (stateVar19)
    temp = [4 19];  
    for i =1:size(temp,2)
      figure, hold on
      plot(tt,squeeze(stateVarFCAll(:,temp(i),:)));
      if meanStateVarFC == 0
      else
        plot(tt,meanStateVarFC(:,temp(i)),'k','LineWidth',2);
      end
      axis([0 65 -inf inf])
      xlabel('Time (min)')
      ylabel(statevar_legend(temp(i)))
    end
end

if plotOutput == 1
    plot_rows = ceil(size(plotOutput_Choice,2)/plot_cols);
    figure('Name','Outputs timeseries (abs)')
    for j = 1:size(plotOutput_Choice,2)
      plotIndex = plotOutput_Choice(j);
      subplot(plot_rows,plot_cols,j), hold on
        plot(tt,squeeze(OutputAll(:,plotIndex,:))); 
      if meanOutput == 0
      else
        plot(tt,meanOutput(:,plotIndex),'k','LineWidth',2);
      end
      title(output_legend(plotIndex))
      if j == 2 || j == 3 || j ==4
        axis([-inf inf 0 1e-3])
      end
    end
       
    % To plot dashed line at 0 in J_F1 graphs (to see where it goes into reverse)
    temp = zeros(1,size(tt(4:stepsize*100:end),2));     
    
    %Plot OCR (Output4), sampled at 6 minute intervals (stepsize = 0.1)
    figure('Name','OCR, JATPcytoProd, J_F1 - raw and converted')
    subplot(3,3,1), hold on
    for i = 1:numSims, plot(tt(4:stepsize*600:end),OutputAll(4:stepsize*600:end,4,i)); end
    if meanOutputFC == 0
    else  % Plot mean trace if multiple simulations were run    
    plot(tt(4:stepsize*600:end),meanOutput(4:stepsize*600:end,4),'k','LineWidth',2);
    end
    axis([0 80 0 10e-4])
    xlabel('Time (min)')
    ylabel('J_C4: OCR (mol O_2/s/l of mito)')
    
    % Plot cytosolic ATP production (Output 27)
    subplot(3,3,4), hold on
    for i = 1:numSims, plot(tt(4:stepsize*100:end),OutputAll(4:stepsize*100:end,27,i)); end
    if meanOutputFC == 0
    else
    plot(tt(4:stepsize*100:end),meanOutput(4:stepsize*100:end,27),'k','LineWidth',2);
    end
    xlabel('Time (min)')
    ylabel('Output27 (J_ATPcytoProd): mol ATP/s/l of mito')
    
    % Plot J_F1 (mitochondrial ATP production?)
    subplot(3,3,7), hold on
    for i = 1:numSims, plot(tt(4:stepsize*100:end),OutputAll(4:stepsize*100:end,5,i)); end
    % Plot line at J_F1 = 0 (to determine when ATP synthase goes in
    % reverse)
    plot(tt(4:stepsize*100:end),temp,'k--')
    if meanOutputFC == 0
    else
    plot(tt(4:stepsize*100:end),meanOutput(4:stepsize*100:end,5),'k','LineWidth',2);
    end
    xlabel('Time (min)')
    ylabel('J_F1: mol ATP/s/l of mito')
           
    % Convert simulations to mol X/min/ug protein to align with experimental data
    dim = 1;
    [OCR_converted,Output27Conv,Output5Conv,...
    J_ATP_total,prop_cyto,prop_mito,convert] = unitsConvert(OutputAll,dim);

    subplot(3,3,2), hold on
    for i = 1:numSims
        Output4Conv = OutputAll(:,4,i)*convert;
        plot(tt(4:stepsize*600:end),Output4Conv(4:stepsize*600:end)); 
    end
    if meanOutputFC == 0
    else
        meanOutput4Conv = meanOutput(:,4)*convert;
        plot(tt(4:stepsize*600:end),meanOutput4Conv(4:stepsize*600:end),'k','LineWidth',2);
    end
    axis([0 80 0 20e-12])
    xlabel('Time (min)')
    ylabel('J_C4: OCR (mol O_2/min/ug protein)')
    
    subplot(3,3,5), hold on
    for i = 1:numSims
        Output27Conv = OutputAll(:,27,i)*convert;
        plot(tt(4:stepsize*100:end),Output27Conv(4:stepsize*100:end));
    end
    if meanOutputFC == 0
    else
        meanOutput27Conv = meanOutput(:,27)*convert;
        plot(tt(4:stepsize*100:end),meanOutput27Conv(4:stepsize*100:end),'k','LineWidth',2);
    end
    xlabel('Time (min)')
    ylabel('Output27 (J_ATPcytoprod): mol ATP/min/ug protein')
    
    subplot(3,3,8), hold on
    for i = 1:numSims
        Output5Conv = OutputAll(:,5,i)*convert;
        plot(tt(4:stepsize*100:end),Output5Conv(4:stepsize*100:end));
    end
    plot(tt(4:stepsize*100:end),temp,'k--')
    if meanOutputFC == 0
    else
        meanOutput5Conv = meanOutput(:,5)*convert;
        plot(tt(4:stepsize*100:end),meanOutput5Conv(4:stepsize*100:end),'k','LineWidth',2);
    end
    xlabel('Time (min)')
    ylabel('J_F1: mol ATP/min/ug protein')
    
    % Plot proportion of cyto ATP production
    subplot(3,3,6), hold on
    for i = 1:numSims
        plot(tt(4:stepsize*100:end), prop_cyto(4:stepsize*100:end,:,i));
    end
    if meanOutputFC == 0
    else
        temp = nanmean(prop_cyto,3);
        plot(tt(4:stepsize*100:end), temp((4:stepsize*100:end)),'k','LineWidth',2);
    end
    axis([-inf inf 0 100])
    xlabel('Time (min)')
    ylabel('Cytosolic ATP production (%)')
    
end

if plotOutputFC == 1
    plot_rows = ceil(size(plotOutputFC_Choice,2)/plot_cols);
    figure('Name','Outputs timeseries (foldchange)')
    for j = 1:size(plotOutputFC_Choice,2)
      plotIndex = plotOutputFC_Choice(j);
      subplot(plot_rows,plot_cols,j), hold on
      for i = 1:numSims
        plot(tt,OutputFCAll(:,plotIndex,i)); 
      end
      if meanOutputFC == 0
      else
        plot(tt,meanOutputFC(:,plotIndex),'k','LineWidth',2);
      end
      title(output_legend(plotIndex))
      %axis([-inf plot_time -inf inf])
    end
    %Plot OCR, sampled at 6 minute intervals (stepsize = 0.1)
    figure, hold on
    for i = 1:numSims, plot(tt(4:stepsize*600:end),OutputFCAll(4:stepsize*600:end,4,i)); end
    if meanOutputFC == 0
    else
    plot(tt(4:stepsize*600:end),meanOutputFC(4:stepsize*600:end,4),'k','LineWidth',2);
    end
    axis([0 65 0 3])
    xlabel('Time (min)')
    ylabel('OCR')
end
end