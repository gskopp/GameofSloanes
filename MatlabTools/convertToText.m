% Convert the Matlab .mat files named dxn_init.mat containing a dxn matrix
% with columns the vectors with variable name F to a .txt file of the required format. 
% Put the files in a subdirectory called tempFileConvert.

packDir=dir('tempFileConvert/*.mat');
for ii = 1:length(packDir)
    packName = packDir(ii).name;
    load(strcat('tempFileConvert/',packName));
    
    [startIndex,endIndex] = regexp(packName,'_[a-zA-Z0-9]+\.');
    
    saveText(F,packName((startIndex+1):(endIndex-1)));
end