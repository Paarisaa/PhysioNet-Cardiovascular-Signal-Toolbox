%	OVERVIEW:
%       This basic demo will allow you to load an ECG file in matlab 
%       compatible wfdb format, detect the locations of the R peaks,
%       perform signal quality (SQI) analysis and plot the results.
%
%   OUTPUT:
%       A figure with the loaded ECG signal and detected peaks will be
%       generated
%   -----------------------------------------------------------------------
%   OUTPUT of the part of the code added by Parisa Sarikhani:
%
%       The variable named "Periodicity" shows the estimed periodicity of
%       the signal in beats/sec.
%        
%       A 3D figure shows the result of lining up the detected peaks of the
%       signal.
%   -----------------------------------------------------------------------
%
%   DEPENDENCIES & LIBRARIES:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   REFERENCE: 
%       Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%       Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian   
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Where are the data, in this demo they are located in a subfolder
InputFolder = [pwd filesep 'TestData' filesep 'mitdb-Arrhythmia']; % path to the folder where you data are located
SigName = '200m';

% load the ecg using rdmat
% [tm,sig,Fs] = rdmat([InputFolder filesep SigName]);
% rdmat was not a defined. I used load instead:
% load the ecg signal using load (it loads a variable called val)
load([InputFolder filesep SigName]);
% the signal has two channels, from now on we will use just one 
ecg = val(1,:);
% Get sampling frequency Fs from header file
sigInfo = readheader2([InputFolder filesep SigName '.hea']);
Fs = sigInfo.freq;
% time vector for visualization (in seconds)
tm = 0:1/Fs:(length(ecg)-1)/Fs;

% plot the signal
figure(1)
plot(tm,ecg);
xlabel('[s]');
ylabel('[mV]')


% Detection of the R-peaks using the jqrs.m function included in the
% toolbox, requires to set initialization parameters calling the
% InitializeHRVparams.m function

HRVparams = InitializeHRVparams('Demo');
% set the exact sampling frequency usign the one from the loaded signal
HRVparams.Fs = Fs;
% call the function that perform peak detection
r_peaks = jqrs(ecg,HRVparams);

% plot the detected r_peaks on the top of the ecg signal
figure(1)
hold on;
plot(r_peaks./Fs, ecg(r_peaks),'o');
legend('ecg signal', 'detected R peaks')



%%
%-------------------------------------------------------------------------%
% The following parts of the code was added by Parisa Sarikhani on
% 09/12/2018 as part of the BMI500 lab assignments.
%   OVERVIEW:
%     This part of the code will estimate the periodicity of the signal in 
%     (beats/sec) and then lines up the detected peaks.
%
%   OUTPUT:
%       The variable named "Periodicity" shows the estimed periodicity of
%       the signal in beats/sec.
%       A 3D figure shows the result of lining up the detected peaks of the
%       signal.
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%% Estimate the periodicity %%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%Calculate the distance of contiguous peaks and convert it from seconds to
%beats per second
diff_peaks=r_peaks(1,2:end)./Fs-r_peaks(1,1:end-1)./Fs;
periodicity=1/(sum(diff_peaks)/length(diff_peaks));%Estimated periodicity(HR(beats/sec)) is simply the inverse of an average over distance of contiguous peaks
display(periodicity)%print the periodicity


%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%% Line Up the Detected Peaks %%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% find the length of left and right intervals (left_length and right_length 
% respectively) of each peak which are half the distance of two contiguous peaks:

for ipeak=1:length(r_peaks) %for each detected peak
    if ipeak==1 %first interval
        left_length(1,ipeak)=r_peaks(1,1);%
        right_length(1,ipeak)=floor((r_peaks(1,2)-r_peaks(1,1))/2);
    elseif ipeak==length(r_peaks) %intermediate intervals
        left_length(1,ipeak)=ceil((r_peaks(1,end)-r_peaks(1,end-1))/2);
        right_length(1,ipeak)=650000-r_peaks(1,end);
    else %last interval
        left_length(1,ipeak)=ceil((r_peaks(1,ipeak)-r_peaks(1,ipeak-1))/2);
        right_length(1,ipeak)=floor((r_peaks(1,ipeak+1)-r_peaks(1,ipeak))/2);
    end
end
%find the maximum length of left and right lengths of the split signals
max_length=max([left_length,right_length]);

%set the location of peaks in each split according to the split with maximum
%length:
peak_location=max_length;

%Making all splits the same size to line up the peaks and create a 3D plot
%showing the result:
figure;
%NaN padding to the beginning and end of each split signal to
%make them of the same size. Each peak is located in the middle of its
%corresponding split

Same_size_split=[NaN*ones(1,max_length-left_length(1,1)),ecg(1:r_peaks(1,1)+floor((r_peaks(1,2)-r_peaks(1,1))/2)+1),...
    NaN*ones(1,max_length-right_length(1,1))];%Nan Padding for the first interval (split ecg signal) 
plot3((1:2*max_length+1)./Fs,1*ones(size(Same_size_split)),Same_size_split)%Plot of the first split

for jpeak=2:length(r_peaks)-1
    hold on
    Same_size_split=[NaN*ones(1,max_length-left_length(1,jpeak)),ecg(r_peaks(1,jpeak)-left_length(1,jpeak):r_peaks(1,jpeak)+right_length(1,jpeak)),...
        NaN*ones(1,max_length-right_length(1,jpeak))];%NaN padding for intermediate intervals
    plot3((1:2*max_length+1)./Fs,jpeak*ones(size(Same_size_split)),Same_size_split)%Plot intermediate split signals
end
hold on
j=2;

Same_size_split=[NaN*ones(1,max_length-left_length(1,end)),ecg(r_peaks(1,end)-left_length(1,end):r_peaks(1,end)+right_length(1,end)),...
    NaN*ones(1,max_length-right_length(1,end))];%NaN padding for the last interval
plot3((1:2*max_length+1)./Fs,(jpeak+1)*ones(size(Same_size_split)),Same_size_split);%plot the last interval
%Add title and axis labels for figure
title('Result of lining up the detected peaks')
xlabel('[S]')
ylabel('Number of split signals')
zlabel('[mV] ')

