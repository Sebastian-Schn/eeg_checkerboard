function [in_dirs, out_dirs, sjs] = prepare_workspace
clc
close all
clear

in_dirs = '';               % path to my_dataset dir 
out_dirs = '';              % path to dir where output will be saved. 
addpath '';                 %add spm to path
spm('defaults', 'eeg');

sjs = {
    'sub-01',...
    'sub-02',...
    'sub-03',...
    'sub-04',...
    'sub-05',...
    'sub-06',...
    'sub-07',...
    'sub-08',...
    'sub-09',...
    'sub-10',...
    'sub-11',...
    };
end
