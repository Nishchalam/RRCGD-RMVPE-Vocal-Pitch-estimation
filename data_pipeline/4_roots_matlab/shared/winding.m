% function w = winding(s, flag)
%
% PURPOSE
% Count the number of times the complex polygonal curve s winds around 0.
%
% REQUIRED INPUT
% s - Any polygonal curve in the complex plane.
% flag = 1 means close the polygon by repeating the last value at the beginning.
%      = 0 means do not close. This is useful if you want to use a real FFT
%        to only evaluate at upper half-plane unit-circle roots instead of using
%        a slower complex FFT.

% Do not use this. Instead compile the mex file in winding.c
%
% AUTHOR
%  James W. Fox

function w=winding(s,flag)
w=0;
