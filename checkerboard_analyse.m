cd '' % path to directory with prepare_workspace.m script
[in_dirs, out_dirs, sjs] = prepare_workspace;

for indx = 1:length(sjs)
    
SJ = sjs{indx};
in_dir = fullfile(in_dirs, sjs{indx}, 'ses-out');
out_dir = fullfile(out_dirs, sjs{indx});

if not(exist(out_dir, 'dir'));
  mkdir(out_dir);
end

S      = [];
S.sdir = in_dir;
S.SJ   = SJ;
S.dataset = fullfile(in_dir, [SJ  '_task-checker_eeg.eeg']);
S.outfile = fullfile(out_dir, [SJ '_task-checker_eeg']);
S.mode                  = 'continuous'; % convert data as continuous 
S.checkboundary         = 1; % 1 = check if there are breaks in the file and do not read across those breaks [default], 0 = ignore breaks (not recommended).
S.saveorigheader        = 0; % 1 = save original data header with the dataset, 0 = do not keep the original header [default]
spm_eeg_convert(S); % Convert various M/EEG formats to SPM12 format

%--------------------------------------------------------------------------
% prepare montaging
S                       = []                            ; % spm_eeg_montage structure initilization
S.D                     = fullfile(out_dir, [SJ '_task-checker_eeg.mat']); % source filename


D                       = spm_eeg_load(S.D)             ; % load MEEG object to obtain channel labels
S.montage.labelorg      = D.chanlabels([1:30 33:64])    ; % N x 1 cell-array: original labels
S.montage.labelnew      = D.chanlabels([1:30 33:64])    ; % N x 1 cell-array: new labels

% create re-reference matrix 
N                       = numel(S.montage.labelorg)     ; % number of channels
ref_tra                 = NaN(N,N)                      ; % reference matrix initialization

% create new montage matrix for  average reference
for i = 1:N
for j = 1:N        
    if i == j 
        ref_tra(i,j) = (N-1)/N;
    else
        ref_tra(i,j) = -1/N;
    end
end
end

S.montage.tra           = ref_tra                       ; % M x N matrix
S.keepothers            =  0                          ; % 'yes'/'no': keep or discard the channels not involved in the montage [default: 'yes']
S.prefix                = 'r_'                          ; % prefix for the output file (default - 'M')

% create new data montage
spm_eeg_montage(S);


% High-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = fullfile(out_dir, ['r_' SJ '_task-checker_eeg.mat']) ; % source filename
S.band                  = 'high'                        ; % filterband [low|high|bandpass|stop]
S.freq                  = 0.1                           ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = 'fh'                         ; % prefix for the output file (default - 'f')
% filter the data
spm_eeg_filter(S);

% Low-Pass Filtering
% -------------------------------------------------------------------------
S.D                     = fullfile(out_dir, ['fhr_' SJ '_task-checker_eeg.mat']); % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 40                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = 'fl'                          ; % prefix for the output file (default - 'f')
%filter the data
spm_eeg_filter(S); 


% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.timewin               = [-100 1000]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = 'e'                               ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

S.D                     = fullfile(out_dir, ['flfhr_' SJ '_task-checker_eeg.mat']);
spm_eeg_epochs(S);

% Artefact correction
% -------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = fullfile(out_dir, ['eflfhr_' SJ '_task-checker_eeg.mat']); % source filename
S.methods.channels      = {'EEG'};
S.methods.fun           = 'peak2peak';
S.methods.settings.threshold = 50;
S.prefix                     = 'a';
spm_eeg_artefact(S);


% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.prefix                = 'm'                           ; % prefix for the output file (default - 'm')

S.D                     = fullfile(out_dir, ['aeflfhr_' SJ '_task-checker_eeg.mat']); % source filename
spm_eeg_average(S);
end

% Grand Average
%-------------------------------------------------------------------------
S = [];
S.D = [fullfile(out_dirs, 'sub-01\maeflfhr_sub-01_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-02\maeflfhr_sub-02_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-03\maeflfhr_sub-03_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-04\maeflfhr_sub-04_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-05\maeflfhr_sub-05_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-06\maeflfhr_sub-06_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-07\maeflfhr_sub-07_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-08\maeflfhr_sub-08_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-09\maeflfhr_sub-09_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-10\maeflfhr_sub-10_task-checker_eeg.mat');...
       fullfile(out_dirs, 'sub-11\maeflfhr_sub-11_task-checker_eeg.mat')];
S.outfile = fullfile(out_dirs, 'gm');
spm_eeg_grandmean(S);

