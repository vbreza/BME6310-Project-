
function [OCR_converted,Output27Conv,Output5Conv,...
    J_ATP_total,prop_cyto,prop_mito,convert] = unitsConvert(Output2Convert,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert simulations to experimental units (pmol X/min/ug protein)
% to align with experimental data for plotting (and calibration)
% dim = dimension of Output2Convert (1 => Output(:,:,:); 2 =>
% Output(:,:,:,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert OCR to experimental units (pmol/min/ug protein)
% Simulation OCR = mol O2 / s / litre of mitochondria
% Experiment OCR = pmol O2 / min / ug protein 

% Simulation J_ATP = mol ATP / s / litre of mitochondria
% Experiment J_ATP = pmol ATP / min / ug protein
    
sec2min = 60;   % to convert from /s to /min

% to convert from /l of mito to total mito volume 
% Ward 2007 CGNs mitochondrial volume = 1.76 +/- 0.15 e-14 litres
litre2mito = 5e-14; 

% to convert from 1 cell to # cells/well 
% 300,000 = Pierre's seeding density
cell2wells = 300000;   

% to convert # cells to ug protein 
% 40-50 ug protein in 300000 cells / well (Pierre's measurements)
cells2protein = 45;     

convert = (sec2min * litre2mito * cell2wells) / cells2protein;

mol2pmol = 1e12;

switch dim
    case 1
    % OCR
    OCR_converted = Output2Convert(:,4,:) * convert * mol2pmol;
    % Cytosolic ATP production
    Output27Conv = Output2Convert(:,27,:) * convert * mol2pmol; 
    % J_F1 converted to mito ATP production
    Output5Conv = Output2Convert(:,5,:) * convert * mol2pmol;   
    
    case 2
    OCR_converted = Output2Convert(:,4,:,:) * convert * mol2pmol;
    Output27Conv = Output2Convert(:,27,:,:) * convert * mol2pmol; 
    Output5Conv = Output2Convert(:,5,:,:) * convert * mol2pmol;   

    case 3
    OCR_converted = Output2Convert(:,4,:,:,:) * convert * mol2pmol;
    Output27Conv = Output2Convert(:,27,:,:,:) * convert * mol2pmol; 
    Output5Conv = Output2Convert(:,5,:,:,:) * convert * mol2pmol;   
end

    % Total ATP production
    J_ATP_total = Output5Conv + Output27Conv;
    % Proportion of total ATP production contributed by cytosolic ATP
    % production
    prop_cyto = (Output27Conv ./ J_ATP_total) * 100;
    % Proportion of total ATP production contributed by mitochondrial ATP
    % production
    prop_mito = (Output5Conv ./ J_ATP_total) * 100;

end