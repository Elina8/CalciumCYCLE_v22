% Created by Evangeline Tzatzalos
% Most recently updated on 05/18/2016
% Copyright (C) of Evangeline Tzatzalos
% For more information contact etzatzalos@gmail.com

function CalciumCYCLE_v22(sampleID, folderName, firstFileToScan, first_sample, smoothingFactor, thresh_maxpeaks, thresh_diff2, MPD, MPD2)
directory='C:/Users/Elina/Desktop/GitHubRepository_1/';
% Script is written for MATLAB(R)
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
%    sampleID                  % name of identification for sample; E### EHM followed by ID, M#### monolayer followed by date.  It will be the first characters in all saved files.
%    folderName                % the name of the ConditionX
%    directory_name            % the highest order folder for the experiment
%    firstFileToScan           % default is 1;     % identifies the first *.csv file to analyze
%    first_sample              % default is 1;     % identifies the first recording in the *.csv file to analyze 
%    smoothingFactor           % default is 6;     % Extent of smoothing, 0 to 12
%    thresh_maxpeaks           % default is 0;     % maximum peaks will be calculated above this value
%    thresh_diff2              % default is 5;     % maximum peaks of 2nd derivatives will be calculated above this value for determining runoff points
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
%    F3. Data files are saved as *.dat under OUTPUT_allRegions.  Each *.dat file contains the following data from a single cell of OUTPUT:
%           Amplitude
%           Time-to-peak
%           Beating Rate
%           TD50
%           TD90
%           Decay tau
%    F4. Data files containing mean and standard deviation are saved as as *.dat under OUTPUT_mean_sd
%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampleID='1';    %E### EHM followed by ID, M#### monolayer followed by date
% folderName = 'ConditionX';
% directory='C:/Users/Elina/Desktop/GitHubRepository_1/';
% firstFileToScan = 1;
% first_sample= 1;
% 
% % Filtering
% smoothingFactor=4;       % removed smoothing and replaced with IIR filter
% 
% % Input Parameters
% thresh_maxpeaks=0;       % Maximum peaks threshold; see below "Amplitude", 0.9
% thresh_diff2=5 ;         % threshhold for maximum 2nd derivatives of intensity curves (FLUOsmooth), 3.5
% MPD=8;                   % minimum peak difference to identify maximums; increase MPD if too many maxs; decrease MPD if too few maxes
% MPD2=5;                  % for finding peaks of the 2nd derivatives for identifying runoff points

% Plotting
ylim_max=5;
ylim_min=0;
% ylim_max_tau=2;
xlim_min=0;
xlim_max=15;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE CELL WITH NORMALIZED DATA FILES (WITHIN CATEGORY, i.e. BaseSpon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(strcat(directory,folderName));
csvFiles = dir('*.csv');     % Identify the CSV files in the folder of ConditionX
numFiles = length(csvFiles); % Count the # of files that are csv in the folder of ConditionX

% Read in all data that fit into the same category (about 4-6 files per category). 
for aa=1:numFiles 
    rawDataCell{aa} = csvread(csvFiles(aa).name,2,3); % Read in csvFiles, read in starting from 2nd row (no header) and 3rd column (time columns)
    TIME{aa}= csvread(csvFiles(aa).name,2,0);         % Cell containing all 
end
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTRODUCING RAW DATA:
% CREATE A CELL WHICH CONTAINS THE MATRIX OF NORMALIZED DATA FOR EACH VIDEO FILE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for bb=firstFileToScan:numFiles; % Each iteration is an EHM area(one video file).  There should be 1-5 EHM areas/video viles.
    % Keep every 6th column of the raw data.  These are the columns with the normalized data
    NormCols=(size(rawDataCell{bb},2))/6;            % Number of normalized columsn.  Also the number of regions within a recorded region/video file.
    normData{bb}=rawDataCell{bb}(:,(1:NormCols)*6);  % Redefine the matrices in the cell to only contain every 6th column (which corresponds to the normalized data) only column 6,12,18...60)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FLUO4 INTENSITY VECTOR (FOR 25 seconds at 25 frames per second (fps))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    numSamples=size(normData{bb},2);
    last_sample=numSamples;
    FLUO=normData{bb}; % FLUO is matrix.  Redefine the Normalized Data Matrix (containining Fluo4 Intensities) to match the below script which was written first.FLUO is a matrix with all normalized intensity values.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE TIME VECTOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nrow, ncol]=size(FLUO);
    frames=TIME{bb}(end,1); 
    totalTime=TIME{bb}(end,2); % last time value in raw data file
    frame_rate=frames/totalTime;
    time=linspace(0,totalTime,frames); %Create a time vector from 1 to totalTime seconds in increments of number of frames
    time = time.'; % make row into column
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF SIGNAL METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=first_sample:last_sample  % Each iteration is a region of the recorded area (video file).  
        samplePosition=i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. AMPLITUDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        FLUOnoNAN=FLUO(~isnan(FLUO(:,i))); %remove NaN values
        minFLUO=min(FLUO(:,samplePosition));
        FLUO(:,samplePosition)=FLUO(:,samplePosition)-minFLUO; %zeroed
        FLUOsmooth=smooth(time,FLUO(1:length(time),samplePosition),smoothingFactor);
        FLUOsmooth=FLUOsmooth(1:length(time));
        FLUOsmoothfloor=FLUOsmooth-min(FLUOsmooth);
        [MAXpeaks,MAXlocs]=findpeaks(FLUOsmooth, 'minPeakDistance', MPD);
        MAXpeaks2=MAXpeaks(MAXpeaks>thresh_maxpeaks);
        MAXlocs2=MAXlocs(MAXpeaks>thresh_maxpeaks);

        %Identify Runoff points (beginning of the episode)
        diff=(FLUOsmooth(2:end)-FLUOsmooth(1:end-1))./(time(2:end)-time(1:end-1));
        diff2=(diff(2:end)-diff(1:end-1))./(time(2:end-1)-time(1:end-2));
        [Diff2peaks,Diff2locs]=findpeaks(diff2,'minPeakDistance', MPD2);

        figure; hold on;
        plot(time,FLUO(1:length(time),samplePosition),'k'); % time is located in 4th column
        plot(time,FLUOsmooth,'r');
%         plot(time(1:end-1),diff,'m-')
%         plot(time(1:end-2),diff2,'c-')
        ylim([ylim_min ylim_max]);
        xlim([xlim_min xlim_max]);
        title(sprintf('Column # = %d of %d',i,ncol));
        legend('Normalized Signal','Filtered Signal')

        Runofflocs = Diff2locs(Diff2peaks>thresh_diff2);   %filter the location of the maximum second derivatives
        Runoffvals = Diff2peaks(Diff2peaks>thresh_diff2);   %filter the values of the maximum second derivatives

        % Pairing MAX and RUNOFF POINTS
        MAXlocs2b=MAXlocs2(MAXlocs2>Runofflocs(1)); % only want the maximum that are after the first runoff point
        Runofflocs2=0;
        for k=1:length(MAXlocs2b)
            for j=1:length(Runofflocs)
                if MAXlocs2b(k)>Runofflocs(j)  % isolate MINpeaks that are to the left of MAXpeaks
                   Runofflocs2(k)=Runofflocs(j);
                end
            end
            if   k==j                              % end if there are more maxes than mins
                 break
            end
        end
        Runofflocs2=Runofflocs2';
        MAXpeaks3=MAXpeaks2(1:k);
        MAXlocs3=MAXlocs2b(1:k);
        % MAXlocs3 and Runofflocs2 are paired
        amplitude=MAXpeaks3(:)-FLUOsmooth(Runofflocs2);
        plot(time(MAXlocs3),MAXpeaks3, 'go'); %FLUO(MAXlocs3,4)
        plot(time(Runofflocs2),FLUOsmooth(Runofflocs2),'bo')
        OUTPUT{bb}(i,1)=mean(amplitude); %find mean of non-zero values
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. TIME-TO-PEAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        timeToPeak=time(MAXlocs3)-time(Runofflocs2);
        OUTPUT{bb}(i,2)=mean(timeToPeak);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. BEATING RATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        time2=time(~isnan(time));  %remove zeros from time vector
        beatingRate=length(MAXpeaks3)/time2(MAXlocs3(end));
        OUTPUT{bb}(i,3)=beatingRate*60;   %beats per minute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. TD50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        point50up_val=zeros(length(MAXpeaks3)-1); % temporary and output matrices should reset for each region
        point50up_time=zeros(length(MAXpeaks3)-1);
        temp50_locs=zeros(length(MAXpeaks3)-1);
        point50down_locs=zeros(length(MAXpeaks3)-1);
        point50down_time=zeros(length(MAXpeaks3)-1);
               
        for m=1:length(MAXpeaks3)-1 % want to look at all but the last peak because the last peak may be incomplete
            [f,p]=fit(time(Runofflocs2(m):MAXlocs3(m)),FLUOsmooth(Runofflocs2(m):MAXlocs3(m)),'linearinterp');    %fit a linear equation, f(x), between the runoff point and the max

            % upswing
            point50up_val(m) = (FLUOsmooth(MAXlocs3(m))+FLUOsmooth(Runofflocs2(m)))/2 ;  % point halfway between the runoff and the maximum
            point50up_dist = @(x) (point50up_val(m) - f(x)); % Goal is to calculate x at 0.5 y (0.5fmax).  Define a function, halfpoint_dist , which compares the halfway point to the fitted function.  
            point50up_time(m) = fzero(point50up_dist, [time(Runofflocs2(m)) time(MAXlocs3(m))]);  % Go through a range of x values (as defined in 2nd parameter), until the the halfpoint_dist equals the fitted function value.  Report that x-value.  fzero is a function that does that comparison.
            % downswing
            locs50 = MAXlocs3(m):length(FLUOsmooth); %temporary locations from the maximum intensity point to the end of the intensity curve
            temp50 = time(locs50(FLUOsmooth(MAXlocs3(m):end)<=point50up_val(m))); %
            point50down_time(m) = temp50(1);
            temp50_locs=locs50(FLUOsmooth(MAXlocs3(m):end)<=point50up_val(m));
            point50down_locs(m)=temp50_locs(1);
        end

        TD50=mean(point50down_time(1:min(size([point50up_time(:,1) point50down_time(:,1)],1)))-point50up_time(1:min(size([point50up_time(:,1) point50down_time(:,1)],1)))); % change made on 1/28/2016    
        OUTPUT{bb}(i,4)=TD50;

        tempMinimum1=min([length(MAXpeaks3) length(point50up_time) length(point50up_val)]);
        tempMinimum2=min([length(MAXpeaks3) length(point50down_time) length(point50up_val)]);
        plot(point50up_time(1:tempMinimum1),point50up_val(1:tempMinimum1),'g*')
        plot(point50down_time(1:tempMinimum2),point50up_val(1:tempMinimum2), 'b*')
        legend('Normalized Signal', 'Filtered Signal','MAXpeaks','MINpeaks','Upswing @ 50% & 90%','Downswing @ 50% & 90%')
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %5. TD90
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        point90up_val=0;    % resent vectors so it doesn't carry over values from previous recording
        point90up_time=0;
        point90down_time=0;
        point90down_locs=0;
        for n=1:length(MAXpeaks3)
            [f,p]=fit(time(Runofflocs2(n):MAXlocs3(n)),FLUOsmooth(Runofflocs2(n):MAXlocs3(n)),'linearinterp');    %fit a linear equation, f(x), between the runoff point and the max

            % upswing
            point90up_val(n) = 0.1*(FLUOsmooth(MAXlocs3(n))-FLUOsmooth(Runofflocs2(n))) + FLUOsmooth(Runofflocs2(n)); % intensity value from beginning of the episode + 10% of episode
            point90up_dist = @(x) point90up_val(n) - f(x); % Goal is to calculate x at 0.1 y (0.1fmax).  Define a function, halfpoint_dist , which compares the halfway point to the fitted function.  
            point90up_time(n) = fzero(point90up_dist, [time(Runofflocs2(n)) time(MAXlocs3(n))]);  % Go through a range of x values (as defined in 2nd parameter), until the the halfpoint_dist equals the fitted function value.  Report that x-value.  fzero is a function that does that comparison.

            locs90 = MAXlocs3(n):length(FLUOsmooth); %temporary locations from the maximum intensity point to the end of the intensity curve
            temp90 = time(locs90(FLUOsmooth(MAXlocs3(n):end)<=point90up_val(n))); %the time at the first point where point falls below 90% point.  if the data drifts, temp90 will be zero.
            if any(temp90)==0        % if temp90 is 0, then point90down_time won't be calculated.  it is better to cut the data at that point.
                break
            else
                point90down_time(n)= temp90(1);
                temp90_locs=locs90(FLUOsmooth(MAXlocs3(n):end)<=point90up_val(n));
                point90down_locs(n)=temp90_locs(1);
            end
        end
    
        point90up_time3=point90up_time(1);   % defining first value
        point90down_time3=point90down_time(1); % defining first value
        for nn=2:min(length(point90down_time), length(point90up_time))-1  % redefine point90down_time to eliminate nonmatched points
            if point90down_time(nn)<point90up_time(nn+1) & point90down_time(nn)>point90down_time(nn-1)
                point90down_time3=[point90down_time3 point90down_time(nn)];
                point90up_time3=[point90up_time3 point90up_time(nn)];
            end
        end
        TD90=mean(point90down_time3-point90up_time3);  % change made on 1/28/2016            

        tempMinimum3=min([length(MAXpeaks3) length(point90up_time) length(point90up_val)]);
        tempMinimum4=min([length(MAXpeaks3) length(point90down_time) length(point90up_val)]);
        plot(point90up_time(1:tempMinimum3),point90up_val(1:tempMinimum3),'g*')
        plot(point90down_time(1:tempMinimum4),point90up_val(1:tempMinimum4), 'b*')

        OUTPUT{bb}(i,5)=TD90;
%         folderName_saveFigs='OUTPUT_analyzed_signal_graphs';    % print figures
%         cd(strcat(directory,folderName_saveFigs));
%         print(strcat(sampleID,'_',folderName,'_',num2str(bb),'-',num2str(i)),'-dpng')
%         cd(strcat(directory,folderName));  % return to previous directory

        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %6. DECAY TAU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         locs=1:length(FLUOsmooth); %location values for FLUOsmooth
         for q=1:min([length(MAXpeaks3) length(point90down_time) length(point90up_time)])-1        %length(MAXpeaks3)-1
            locs_temp=locs(MAXlocs3(q):(end-1)); %locations from peak to end
            temp_EpisodeEnd_locs=locs_temp(diff(MAXlocs3(q):locs(end-1))>=0); %where 1st differential is greater or equal than 0
            pointEpisodeEnd_loc(q)=temp_EpisodeEnd_locs(1); %End of episode locations (first point where 1st differential is 0).
            pointEpisodeEnd_time(q)=time(pointEpisodeEnd_loc(q)); % End of episode time

            for r=MAXlocs3(q):Runofflocs2(q+1)  %Calculate the beginning of the decay region
                if time(r)>=point90down_time(q)
                break
                end
            end
            if point50down_locs(q)>=pointEpisodeEnd_loc(q)      % a missing Runoff point may make point50down_locs larger than pointEpisodeEnd-Loc.  then code fails when it tries to fit this range
                break
            end

            [f_exp,p]=fit(time(point50down_locs(q):pointEpisodeEnd_loc(q)),FLUOsmoothfloor(point50down_locs(q):pointEpisodeEnd_loc(q)),'exp1');
            c=coeffvalues(f_exp);
            tau(q)=-1/c(2);
        end

%         figure;hold on;
%         plot(time,FLUOsmoothfloor(1:length(time)),'Color',[0 0 0]);
%         plot(time(point90down_locs(q):pointEpisodeEnd_loc(q)),f_exp(time(point90down_locs(q):pointEpisodeEnd_loc(q))), 'r*');
%         ylim([0 ylim_max_tau]);
%         hold off;

        OUTPUT{bb}(i,6)=mean(tau);
        OUTPUT{bb}(i,7:11)=[smoothingFactor; thresh_maxpeaks; thresh_diff2; MPD; MPD2];

        % Save figures into OUTPUT_analyzed_signal_graphs
        folderName_saveFigs='OUTPUT_analyzed_signal_graphs';    % print figures
        cd(strcat(directory,folderName_saveFigs));
        print(strcat(sampleID,'_',folderName,'_',num2str(bb),'-',num2str(i)),'-dpng')
        cd(strcat(directory,folderName));  % return to previous directory
%         
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %6. AVERAGE REGIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTPUT_allRegions=OUTPUT{1}(:,1:6);
for cc=firstFileToScan+1:numFiles
    OUTPUT_allRegions=[OUTPUT_allRegions; OUTPUT{cc}(:,1:6)];
end

OUTPUT_mean=mean(OUTPUT_allRegions);
OUTPUT_stdev=std(OUTPUT_allRegions);
OUTPUT_mean_sd=[OUTPUT_mean; OUTPUT_stdev];

OUTPUT_input_parameters=OUTPUT{1}(1,7:11);    % all input parameters (same for each row so you don't have to average them)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %6. SAVE OUTPUT MATRIX
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(strcat(directory, 'OUTPUT_allRegions'));
filename=strcat(sampleID,'_',folderName(1:6),'_All','.dat');
csvwrite(filename,OUTPUT_allRegions)

cd(strcat(directory, 'OUTPUT_mean_sd'));
filename=strcat(sampleID,'_',folderName(1:8), '.dat');
csvwrite(filename,OUTPUT_mean_sd)

cd(strcat(directory, 'OUTPUT_input_parameters'));
filenameInputParameters=strcat(sampleID,'_',folderName(1:6),'_input_parameters','.dat');
csvwrite(filenameInputParameters,OUTPUT_input_parameters)

cd(strcat(directory,folderName));



%% Troubleshooting
%Max peaks (green circle):  change the MPD
%Min peaks (blue circle): change maxdiff2 and MPD2
%When not all maxs and mins (right half of the trace), reduce MPD2
% Remove recordings that display the following features:  
%    a. non-contiguous
%    b. drifting
%    c. interupted transients
%    d. <0.1 amplitude
%    e. when paced beating rates are not near 60bpm (for 1 Hz pacing)

end