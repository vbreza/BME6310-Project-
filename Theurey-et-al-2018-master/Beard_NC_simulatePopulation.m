function Beard_NC_simulatePopulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code adapted from crosst.m of Huber 2011 (PMID: 21364572) 
%%% and regulation.m of Huber 2012 (PMID: 22218564). 
%%% Model originally developed in Beard 2005 (PMID: 16163394)
%%% This code simulates a population for one expt type
%%% - Runs multiple times, varying parameters by +/- X% each time
%%% Can simulate all drug additions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define model parameters
xpar = define_model_parameters;
otherpar = define_other_parameters;

%%% Define simulation time-frames and other time settings
printRequest = 0;   %Set to 1 if you want to print info on time settings
[t_prior, t_start,t_final,t_no_time,stepsize,tt,time]...
    = defineSimulationTimeFrames(printRequest);

protons_per_ATPase = 3;     % three protons per ATP synthase (for output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set what figures to plot
plot_ss = 0;    % Set to 1 to plot steady-state simulations
plotStateVar = 1;  % Set to 1 to plot all raw state variables
plotStateVar_Choice = [4 8 13 19 23 24];   % Choose which outputs you want to plot (listed in getLegends.m)
plotStateVarFC = 0;  %...state variable fold changes
plotStateVarFC_Choice = [4 8 13 19 23 24]; 
plotOutput = 1;       % ...raw outputs  
plotOutput_Choice = [1 2 3 4 5 6 13 15 17 27];   
plotOutputFC = 0;  %...output fold changes
plotOutputFC_Choice = [1 2 3 4 5 6 13 15 17 27]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctot0       = xpar(2);           % Cyt-c from Beard. Total IMS Cyt-c, Cred+Cox, molar
Qtot0       = xpar(3);           % total IMS Ubiquinol, Q+QH2, molar
ADTP_tot    = xpar(4);           % total Adenosine phosphates in cytosol, calculated from Evans_77

%%% Basal Respiratory State
state_fact  = 10/11; 
ATP_e       = state_fact*ADTP_tot;  % ATP level
ADP_e       = ADTP_tot-ATP_e;       % ADP level
fprintf('Respiratory State (ATP:ADP): %0.1i\n', state_fact)
fprintf('Glycolytic Capacity (K_ADTP_dyn): %0.1i\n', otherpar(4))

xo_single_cell  = initial(ADP_e, ATP_e, Ctot0, Qtot0); 

%%% ODE options
options = odeset('RelTol',1e-5, 'AbsTol',1e-8, 'MaxStep',10e-1, ...
    'InitialStep',1e-1, 'MaxOrder',5, 'BDF','on');

%%% Assign legends for output graphs 
[statevar_legend, output_legend] = getLegends; 

%%% Set seed of random number generator (for reproducible simulations)
setRNG = 1;
fprintf('\n********************')
fprintf('\n********************')
if setRNG == 1
    rng(1)
    fprintf('\nInitialised random number generator seed [rng(1)]')
end

%%% Set number of simulations (= number of cells in population)
numSims = 10;
printSim = 1;   % Set to 1 to output simulation numbers to screen
fprintf('\n%i simulations', numSims)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set conditions to simulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define experiment to simulate (based on input, expt)
% 1: Rotenone, Oligo
% 2: AA, Oligo
% 3: Seahorse (Oligo, FCCP, Rot+AA)
% 4: Increased energy demand (IncEnDem, Oligo, FCCP, Rot+AA) 
% 5: FCCP, Rotenone 
expt = 3;
printRequest = 1; %Set to 1 if you want to print expt info
fprintf('\nRunning expt %i', expt)
[rotenone,AA,oligo,CIV,FCCP,energy] = defineExptsToSimulate(expt,t_no_time,printRequest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set impairments to simulate
% Set parameters to 1 for 100% condition, and to <1 to
% simulate impairment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_xpar7 = 1;%0.5;%0.42;%    % x_DH impairment 
v_xpar9 = 1;%2e-1;%9e-2;%4.5e-2;%      % x_C1  
v_xpar10 = 1;%1.2e-1;%7.2e-2;%5.2e-2;%     % x_C3
v_xpar13 = 1;%6e-3;%3.2e-3%2.2e-3;%     % x_C4 impairment
v_xpar15 = 1;%1e-4;%6e-5%4.2e-5;%    % x_F1 impairment
v_xpar21 = 1;%1.5;%2.0;%     % x_Hle  
v_otherpar4 = 1;%0.3;%0.82;%  % K_ADTP_dyn impairment
v_otherpar6 = 1; % K_ADTP_cons (redundant)

if v_xpar7 ~= 1, fprintf('Impairing x_DH [xpar(7)=%0.2i]\n',v_xpar7), end
if v_xpar9 ~= 1, fprintf('Impairing x_C1 [xpar(9)=%0.2i]\n',v_xpar9), end
if v_xpar10 ~= 1, fprintf('Impairing x_C3 [xpar(10)=%0.2i]\n',v_xpar10), end
if v_xpar13 ~= 1, fprintf('Impairing x_C4 [xpar(13)=%0.2i]\n',v_xpar13), end
if v_xpar15 ~= 1, fprintf('Impairing x_F1 [xpar(15)=%0.2i]\n',v_xpar15), end
if v_xpar21 ~= 1, fprintf('Impairing x_Hle [xpar(21)=%0.2i]\n',v_xpar21), end
if v_otherpar4 ~= 1, fprintf('Impairing K_ADTP_dyn [otherpar(4)=%0.2i]\n',v_otherpar4), end

% Prompt user to cancel if incorrect parameters are altered.
fprintf('Have you set breakpoint??\n')
input('All OK? (Ctrl+C if not)\n')
fprintf('********************\n')
fprintf('********************\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run multiple simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xparAll, stateVarAll,stateVarFCAll,stateVarAll1,...
    stateVarAll8,stateVarAll30,stateVarAll45,...
    stateVarAll70, OutputAll,OutputFCAll,OutputAll1,...
    OutputAll8,OutputAll30,OutputAll45,OutputAll70,...
    median_stateVarAll,median_stateVarFCAll,median_stateVarAll1,...
    median_stateVarAll8,median_stateVarAll30,median_stateVarAll45,...
    median_stateVarAll70,median_OutputAll,median_OutputFCAll,...
    median_OutputAll1,median_OutputAll8,median_OutputAll30,...
    median_OutputAll45,median_OutputAll70,...
    mean_stateVarAll,mean_stateVarFCAll,mean_stateVarAll1,...
    mean_stateVarAll8,mean_stateVarAll30,mean_stateVarAll45,...
    mean_stateVarAll70, mean_OutputAll,mean_OutputFCAll,mean_OutputAll1,...
    mean_OutputAll8,mean_OutputAll30,mean_OutputAll45,mean_OutputAll70] ...
    = populationRun(numSims, oligo,rotenone,AA,CIV,FCCP,energy,time,...
    v_xpar7,v_xpar9,v_xpar10,v_xpar13,v_xpar15,v_xpar21,...
    v_otherpar4,v_otherpar6,...
     t_prior,t_start,t_final,tt,stepsize,...
     options, xo_single_cell, plot_ss, statevar_legend, printSim);

% Print again just to make sure you don't miss it!
fprintf('\n********************')
if v_xpar7 ~= 1, fprintf('\nImpairing x_DH [xpar(7)=%0.2i]\n',v_xpar7), end
if v_xpar9 ~= 1, fprintf('\nImpairing x_C1 [xpar(9)=%0.2i]\n',v_xpar9), end
if v_xpar10 ~= 1, fprintf('\nImpairing x_C3 [xpar(10)=%0.2i]\n',v_xpar10), end
if v_xpar13 ~= 1, fprintf('\nImpairing x_C4 [xpar(13)=%0.2i]\n',v_xpar13), end
if v_xpar15 ~= 1, fprintf('\nImpairing x_F1 [xpar(15)=%0.2i]\n',v_xpar15), end
if v_xpar21 ~= 1, fprintf('\nImpairing x_Hle [xpar(21)=%0.2i]\n',v_xpar21), end
if v_otherpar4 ~= 1, fprintf('\nImpairing K_ADTP_dyn [otherpar(4)=%0.2i]\n',v_otherpar4), end
if v_otherpar6 ~= 1, fprintf('\nImpairing K_ADTP_cons [otherpar(6)=%0.2i]\n',v_otherpar6), end
fprintf('********************\n')

if setRNG == 1
    fprintf('\nRandom number generator seed was reset [rng(1)] prior to simulations.\n')
end
fprintf('\nRunning expt %i\n', expt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot timeseries of state variables & outputs
%%% Only plots what you have chosen at beginning of file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotTimeSeries(tt, stepsize, numSims,...
    plotStateVar, plotStateVar_Choice,stateVarAll,...
    plotStateVarFC, plotStateVarFC_Choice, stateVarFCAll,...
    plotOutput, plotOutput_Choice, OutputAll,...
    plotOutputFC, plotOutputFC_Choice, OutputFCAll)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output specific values at specific timepoints to obtain data for 
%%% Theurey et al Fig 2D and 4A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TMRM
tempTMRM = squeeze(OutputAll(:,13,:));    % \Delta\Psi_m
TMRMb = tempTMRM(2/stepsize,:)';    % Baseline
TMRMd1 = tempTMRM(30/stepsize,:)';    % After 1st drug addition
TMRMd2 = tempTMRM(45/stepsize,:)';   % After 2nd drug addition
tempTMRM = [TMRMb TMRMd1 TMRMd2];

tempFCTMRM = squeeze(OutputFCAll(:,13,:));  % Foldchange
TMRMFCd1 = tempFCTMRM(30/stepsize,:)';
TMRMFCd2 = tempFCTMRM(45/stepsize,:)';
tempFCTMRM = [TMRMFCd1 TMRMFCd2];

% NADH
tempNADH = squeeze(stateVarAll(:,4,:)); % NADH
NADHb = tempNADH(2/stepsize,:)';
NADHd1 = tempNADH(30/stepsize,:)';
NADHd2 = tempNADH(45/stepsize,:)';
tempNADH = [NADHb NADHd1 NADHd2];

tempFCNADH = squeeze(stateVarFCAll(:,4,:));
NADHFCd1 = tempFCNADH(30/stepsize,:)';
NADHFCd2 = tempFCNADH(45/stepsize,:)';
tempFCNADH = [NADHFCd1 NADHFCd2];

% OCR
tempOCR = squeeze(OutputAll(:,4,:));   % OCR (raw units)
OCRb = tempOCR(2/stepsize,:)';
OCRd1 = tempOCR(30/stepsize,:)';
OCRd2 = tempOCR(45/stepsize,:)';
tempOCR = [OCRb OCRd1 OCRd2];

tempFCOCR  = squeeze(OutputFCAll(:,4,:));
OCRFCd1 = tempFCOCR(30/stepsize,:)';
OCRFCd2 = tempFCOCR(45/stepsize,:)';
tempFCOCR = [OCRFCd1 OCRFCd2];

% Convert OCR to experimental units
% convert units first
Output2Convert = OutputAll;
dim = 1;
[OCR_converted,Output27Conv,Output5Conv,...
J_ATP_total,prop_cyto,prop_mito,convert] = unitsConvert(Output2Convert,dim);

temp_OCRconv = squeeze(OCR_converted); %(Want timeseries for OCR data)
basal_conv = temp_OCRconv(2/stepsize,:)';    
drug1_conv = temp_OCRconv(30/stepsize,:)';    
drug2_conv = temp_OCRconv(45/stepsize,:)';   
temp_OCRconv = [basal_conv drug1_conv drug2_conv];

% J_Hle
tempJHle = squeeze(OutputAll(:,15,:));
JHleb = tempJHle(2/stepsize,:)';
JHled1 = tempJHle(30/stepsize,:)';
JHled2 = tempJHle(45/stepsize,:)';
tempJHle = [JHleb JHled1 JHled2];

tempFCJHle = squeeze(OutputFCAll(:,15,:));
JHleFCd1 = tempFCJHle(30/stepsize,:)';
JHleFCd2 = tempFCJHle(45/stepsize,:)';
tempFCJHle = [JHleFCd1 JHleFCd2];

% J_F1
tempJF1 = squeeze(OutputAll(:,5,:));
JF1b = tempJF1(2/stepsize,:)';
JF1d1 = tempJF1(30/stepsize,:)';
JF1d2 = tempJF1(45/stepsize,:)';
tempJF1 = [JF1b JF1d1 JF1d2];

tempFCJF1 = squeeze(OutputFCAll(:,5,:));
JF1FCd1 = tempFCJF1(30/stepsize,:)';
JF1FCd2 = tempFCJF1(45/stepsize,:)';
tempFCJF1 = [JF1FCd1 JF1FCd2];

% Set breakpoint here if you want to extract specific information
end