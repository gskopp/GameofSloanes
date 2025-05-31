% Convert the txt files to mat files with synthesis matrices of the packings.  
% Put the files in a subdirectory called tempFileConvert.

packDir=dir('tempFileConvert/*.txt');
for ii = 1:length(packDir)
    packName = packDir(ii).name;
    fid = fopen(strcat('tempFileConvert/',packName),'r');
    
    [startLabelIndex,endLabelIndex] = regexp(packName,'_[a-zA-Z0-9]+\.');
    [~,enddIndex] = regexp(packName,'\d+x');
    d = str2num(packName(1:(enddIndex-1)));
    [startnIndex,endnIndex] = regexp(packName,'x\d+');
    n = str2num(packName((startnIndex+1):endnIndex));

    F=fscanf(fid,'%f');

    F=reshape(F(1:length(F)/2)+1i*F(length(F)/2+1:end),d,n);

    F=F*diag(1./sqrt(diag(F'*F)));

    fclose(fid);
    
    saveF(F,packName((startLabelIndex+1):(endLabelIndex-1)));
end