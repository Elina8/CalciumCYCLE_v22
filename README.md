# CalciumCYCLE_v22
Calculation of quantitative parameters of signals from periodic calcium transients

% Created by Evangeline Tzatzalos

% Most recently updated on 05/18/2016

% Copyright (C) of Evangeline Tzatzalos

% For more information contact etzatzalos@gmail.com

% A. GOAL: to calculate the following from time-dependent fluorescent activity:
%    A1. Amplitude
%    A2. Time-to-peak
%    A3. Beating Rate
%    A4. TD50
%    A5. TD90
%    A6. Decay tau
%
% B. DATA STRUCTURE
%    B1. The highest order folder is labeled according to experiment (directory_name)
%    B2. The following are required embedded folders
%           OUTPUT_allRegions
%           OUTPUT_analyzed_signal_graphs
%           OUTPUT_input_parameters
%           OUTPUT_mean_sd
%           ConditionX              % A new folder should be created for all conditions. 
%    B3. All *.csv files should be placed inside the appropriate ConditionX folder.
%
% C. INPUTS: DATA FILES
%    C1. Input data files are created by exported video files from Axiovision vxxx.
%    C2. Input data files are *.csv. 
%    C3. Time vector is in Column B and must be in "number" format.
%    C4. Normalized fluorescence data of the first sample is in Column I
%    C5. Normalized fluroescence data of every sample thereafter is in every 6th column thereafter (i.e. Column I, O, U, AA...)
%    C6. Column C should be deleted or converted to number format.  It provides information on frame rate which is not necessary.% C. OUTPUT
%
% D. INPUTS: REQUIREMENTS
%    D1. continuous
%    D2. at least 2 peaks
%    D3. at least 10 sec of recording
%    D4. No more than 25 sec of recording
%    D5. no drift from baseline
%    D6. 2 minimums need to flank the selected signal% 
%    D7. Discard samples with arrhythmias or incomplete repolarization of calcium transients
%    D8. Discard samples with high noise
%    D9. Discard with incomplete depolarizations
%
% E. INPUTS: VARIABLES
%    sampleID                  % name of identification for sample; E### EHM followed by ID, M#### monolayer followed by date.  It    %                                will be the first characters in all saved files.
%    folderName                % the name of the ConditionX
%    directory_name            % the highest order folder for the experiment
%    firstFileToScan           % default is 1;     % identifies the first *.csv file to analyze
%    first_sample              % default is 1;     % identifies the first recording in the *.csv file to analyze 
%    smoothingFactor           % default is 6;     % Extent of smoothing, 0 to 12
%    thresh_maxpeaks           % default is 0;     % maximum peaks will be calculated above this value
%    thresh_diff2              % default is 5;     % maximum peaks of 2nd derivatives will be calculated above this value for   
%                                determining runoff points
%    MPD                       % default is 5;     % mean peak distance to identify maximums  
%    MPD2                      % default is 5;     % mean peak distance to identify runoff points
%
% D. OUTPUT
%    D1. The variable OUTPUT is an array of cells, each containing the analytical output from all samples in each .csv file.
%           Column 1: Amplitude
%           Column 2: Time-to-peak
%           Column 3: Beating Rate
%           Column 4: TD50
%           Column 5: TD90
%           Column 6: Decay Tau
%           Column 7: smoothingFactor
%           Column 8: thresh_maxpeaks
%           Column 9: thresh_diff2
%           Column 10: MPD
%           Column 11: MPD2
%
% F. SAVING FEATURES
%    F1. Input parameters are saved as *.dat under OUTPUT_input_parameters
%    F2. Figures are saved as .fig under OUTPUT_analyzed_signal graphs
%           Amplitude
%           Time-to-peak
%           Beating Rate
%           TD50
%           TD90
%           Decay tau
%    F3. Data files are saved as *.dat under OUTPUT_allRegions.  Each *.dat file contains the following data from a single cell of 
%           OUTPUT:
%           Amplitude
%           Time-to-peak
%           Beating Rate
%           TD50
%           TD90
%           Decay tau
%    F4. Data files containing mean and standard deviation are saved as as *.dat under OUTPUT_mean_sd
%%
