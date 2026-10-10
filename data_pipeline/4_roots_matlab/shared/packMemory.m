% function packMemory
%
% PURPOSE:
%	Do garbage collection.
%
% AUTHOR:
%	James W. Fox 1/9/2006

% This is straight from the Matlab documentation for pack
cwd = pwd;
cd(tempdir);
pack
cd(cwd)
clear cwd
