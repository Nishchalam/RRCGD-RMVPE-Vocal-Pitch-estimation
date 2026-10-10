% function rootsIni
%
% PURPOSE:
%   Set, folder, the directory which will contain roots.his and defl.his.
%   This is passed via global variables to deflateFFT.m
%   Set RAM shich fixes MAX_ARRAY_LENGTH_EXPONENT. Both are used in findMultipliers.
%   MAX_ARRAY_LENGTH_EXPONENT is also used in lroots and nroots.
%   DIAGNOSTICS is used in gridsearch, laguerrre3, & unfactor.

% AUTHOR
%  James W. Fox

function rootsIni

% Declare globals
global HISTORYFILE DEFLFILE DIAGNOSTICS RAM MAX_ARRAY_LENGTH_EXPONENT nCPUs;

% Set globals
% ORIGINAL CODE required you to hand-edit this file to point 'folder' at
% an existing writable directory (was hardcoded to 'c:\poly\', Windows-only).
% FIXED: auto-derive a portable 'logs' folder next to this package and
% create it if it doesn't exist yet, so no manual editing is required.
folder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'logs');
if ~exist(folder, 'dir')
    mkdir(folder);
end
folder = [folder, filesep];
HISTORYFILE = [folder, 'roots.his']; % Messages - error or success.
DEFLFILE    = [folder, 'DeflatedRadii.his']; % Radii of deflated roots

DIAGNOSTICS = 0; % Set to 1 to get diagnostic messages in grid search, polishing, and unfactoring.
                  % Useful when debugging a single long factorization; leave 0 for batch runs
                  % over many frames (e.g. the audio benchmark script) to avoid log spam.
RAM = 2; % Physical RAM in GB. Warning! For Windows XP and before, do not set this > 2 MB. Windows only allows 2 GB per process
nCPUs=1; % The number of CPUs on your cluster. The current, 3/3/2006, code requires it to be 1. Hopefully version 3 will be different.

% Compute the maximum length of a real array, 2^MAX_ARRAY_LENGTH_EXPONENT. 
% Note, complex arrays can only be half as long. However, the split radix real rdft produces 
% complex arrays that are 1 + .5*(the length of the real array input).
if(RAM <= .25)
   MAX_ARRAY_LENGTH_EXPONENT = 23; 
elseif(RAM <= .5)
   MAX_ARRAY_LENGTH_EXPONENT = 24; 
elseif(RAM <= 1)
   MAX_ARRAY_LENGTH_EXPONENT = 25; 
elseif(RAM <= 2)
   MAX_ARRAY_LENGTH_EXPONENT = 26;
elseif(RAM <= 4)
   MAX_ARRAY_LENGTH_EXPONENT = 27;
elseif(RAM <= 8)
   MAX_ARRAY_LENGTH_EXPONENT = 28;
elseif(RAM <= 16)
   MAX_ARRAY_LENGTH_EXPONENT = 29;
end

ver = str2double(regexp(version,'^\d+','match','once'));

if(ver < 6) % Matlab 5 can not go as high.
  MAX_ARRAY_LENGTH_EXPONENT = MAX_ARRAY_LENGTH_EXPONENT - 1;
end

% Check to make sure that "folder" is set to an existing folder.
if ~exist(folder, 'dir')
  error('<<Please edit RootsIni.m in folder shared and set the variable "folder" to be a folder on your hard drive.>>');
end

% BATCH MODE: disabled diagnostic output (removed parfor race condition on tempfile deletion)
% tempfile = [folder,'3141592654'];
% fp = fopen(tempfile,'w');
% if(fp == -1)
%   error('<<Please edit RootsIni.m in folder shared and set the variable "folder" to be a folder on your hard drive.>>');
% else
%   fclose(fp);pause(0);
%   delete(tempfile);
% end

tic;
