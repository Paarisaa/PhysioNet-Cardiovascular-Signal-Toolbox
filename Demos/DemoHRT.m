%	OVERVIEW:
%       Demo for Heart Rate Turbulence (HRT) analysis using the Physionet
%       HRV toolbox for Matlab
%       Provided data are from the MIT Physionet Arrhytmia Database 
%       (https://www.physionet.org/physiobank/database/mitdb/)
%
%   OUTPUT:
%       HRT Metrics exported to .cvs files
%
%   REFERENCE: Bauer, Axel, et al. "Heart rate turbulence: standards of 
%              measurement, physiological interpretation, and clinical use:
%              International Society for Holter and Noninvasive 
%              Electrophysiology Consensus." Journal of the American 
%              College of Cardiology 52.17 (2008): 1353-1365.
%	REPO:       
%       https://github.com/cliffordlab/Physionet-HRV-toolbox-for-MATLAB
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian   
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Remove old files generated by this demo
OldFolder = [pwd, filesep, 'OutputData', filesep, 'ResultsHRT'];
if exist(OldFolder, 'dir')
    rmdir(OldFolder, 's');
    fprintf('Old Demo Folder deleted \n');
end


% Initialize settings for demo
HRVparams = InitializeHRVparams('demoHRT'); 
HRVparams.af.on = 0; % No AF analysis for this demo
HRVparams.MSE.on = 0; % No MSE analysis for this demo
HRVparams.DFA.on = 0; % No DFA analysis for this demo
HRVparams.timedomain.on= 0 ;  % No Time Domain analysis for this demo
HRVparams.freq.on= 0 ;  % No Frequency Domain analysis for this demo
HRVparams.Entropy.on = 0 ;  % No Frequency Domain analysis for this demo
HRVparams.poincare.on = 0 ;  % No Frequency Domain analysis for this demo
HRVparams.prsa.on = 0 ;  % No Frequency Domain analysis for this demo
HRVparams.HRT.GraphOn = 1;


% Load qrs and annotations from file
[qrs,ann] = read_ann(['TestData', filesep, 'Physionet_nsr2db' filesep 'nsr004'],'ecg');

% Remove non-rhythm annotations
qrs(ann=='+' | ann=='~' | ann=='|' | ann=='x') = [];
ann(ann=='+' | ann=='~' | ann=='|' | ann=='x') = [];


RRInts = diff(qrs./HRVparams.Fs); % in seconds
tRRInts = qrs(2:end)./HRVparams.Fs;
Labels = ann(2:end);

% Perform HRT analysis
[results, resFilename] = Main_HRV_Analysis(RRInts,tRRInts,'RRIntervals',HRVparams,'nsr004',Labels);

% 3. Compare generated output file with the reference one
currentFile = [HRVparams.writedata filesep resFilename.HRT '.csv'];
referenceFile = ['ReferenceOutput' filesep 'HRTtest.csv'];
testHRV = CompareOutput(currentFile,referenceFile);

if testHRV
    fprintf('** Heart Rate Turbulence Demo: TEST SUCCEEDED ** \n ')
else
    fprintf('** Heart Rate Turbulence Demo: TEST FAILED ** \n')
end

