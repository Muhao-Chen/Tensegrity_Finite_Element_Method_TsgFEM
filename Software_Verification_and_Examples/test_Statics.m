%% test Statics
% Below is to run the statics examples automatically
%%
clear; close all; clc;
%%
testSet='all';
testMode='test';
approveQuestion=0;
startLoc=1;
%% running tests

testFileList={'Main_tower_static.m','Main_folding_3d_Dbar.m','Main_lander.m','main_Jasen_linkage_3.m'};
testFolderList={'01_tower','02_3D_Dbar','03_lander','04_Jasen_linkage\example3'}
%% loop over all examples
for q_test=1:1:numel(testFileList)
    
    fileMessage=['testStatics -> Test file:',num2str(q_test),' of ',num2str(numel(testFileList)),' ',testFileList{q_test}];
    disp(' ');
    disp(fileMessage);
    disp(' ');
    
    % Make testFolder current directory
    testFolder=fullfile(fileparts(mfilename('fullpath')),'Statics_Examples',testFolderList{q_test});
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