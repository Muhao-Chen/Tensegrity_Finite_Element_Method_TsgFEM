%% Automated dynamics test 

% This automated test code was derived from codes from the GIBBON project.
% The source code is here: 
% https://github.com/gibbonCode/GIBBON/blob/master/lib/testGibbon.m
% The license for GIBBON is given below. 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% And the copyright information is:
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON
% contributors.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% Please include their licence and Copyright information when you use this
% file. 

%% test Dynamics
% Below is to run the dynamics examples automatically.
%%
clear; close all; clc;

%% running tests
originFolder=fileparts(mfilename('fullpath'));
testFileList={'main_tower_dynamic_simple.m',...
              'main_10seg_truss_simple.m',...
              'main_double_pendulum_analytical.m',...
              'main_double_pendulum_simple.m',...
              'main_lander.m',...
              'main_tbar_dynamics'...
              };
testFolderList={'01_tower',...
                '02_cantilever_truss',...
                ['03_double_pendulum',filesep,'Analytical_Solution'],...
                ['03_double_pendulum',filesep,'FEM_solution'],...
                '04_Lander',...
                '05_T_bar_Dynamics',...
                };

%% loop over all examples
for q_test=1:1:numel(testFileList)

    fileMessage=['testDynamics -> Test file:',num2str(q_test),' of ',num2str(numel(testFileList)),' ',testFileList{q_test}];
    disp(' ');
    disp(fileMessage);
    disp(' ');

    % Make testFolder current directory
    testFolder=fullfile(fileparts(mfilename('fullpath')),'Dynamics_Examples',testFolderList{q_test});
    addpath(testFolder);
    cd(testFolder);

    mFileNow=fullfile(testFolder,testFileList{q_test});

    initialVars_publish = who;
    save('tempPub.mat'); %Save current variables

    % test the example
    run(mFileNow);

    load('tempPub.mat'); %Load variables
    delete('tempPub.mat'); %Clean up

    choice = questdlg([fileMessage,'. Done. Do you want to proceed?'],testFileList{q_test},'Yes','No','Yes');
    switch choice
        case 'Yes'

        case 'No'
            edit(mFileNow);
            break
    end

    clearvars('-except',initialVars_publish{:});
    close all;
end
cd(originFolder);