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
testFolder=fullfile(fileparts(fileparts(mfilename('fullpath'))),'Statics_Examples');
testFileList={'Main_tower_static.m'};
testFolderList={'01_tower'}
%% loop over all examples
for q_test=1:1:numel(testFileList)
    
    fileMessage=['testStatics -> Test file:',num2str(q_test),' of ',num2str(numel(testFileList)),' ',testFileList{q_test}];
    disp(' ');
    disp(fileMessage);
    disp(' ');
    
    % Make testFolder current directory
    addpmaath(testFolder);
    cd(testFolder); 
    
        mFileNow=fullfile(testFolder,testFileList{q_test});
        
        initialVars_publish = who;
    save('tempPub.mat'); %Save current variables
    
        % test the example
       run(mFileNow);  
       
    load('tempPub.mat'); %Load variables
    delete('tempPub.mat'); %Clean up
    
        clearvars('-except',initialVars_publish{:});
    close all;
end