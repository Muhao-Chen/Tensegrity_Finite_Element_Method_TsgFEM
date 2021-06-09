%% ^_^ Welcome to Tensegrity Finite Element Method(TsgFEM) software! ^_^ %%
% SETUP file to be run only the first time
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

%% Add all necessary functions to MATLAB path%% set up the workspace
clear all;
close all;
clc;
% add the function libraries
addpath( genpath('Function_Library') );
% add the Software Verification and Examples
addpath( genpath('Software_Verification_and_Examples') );
% add User Guide
addpath( genpath('User_Guide') );
% add Videos folder
addpath( genpath('Videos') );
% add JOSS Paper
addpath( genpath('JOSS_Paper') );

%% Open the User_Guide
% cd User_Guide;
% open('User_Guide_Tensegrity_Finite_Element_Method_(TsgFEM).pdf');
% cd ..;
disp('Welcome! Please follow the step-by-step instructions from the User Guide.');