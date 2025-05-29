function F = textToMatlab( fileName, d, n )
%textToMatlab 
%
% Convert a text file, formatted as a long column of first the
% real components of each vector and then the imaginary parts, to a Matlab file.
%
% INPUTS
% fileName: sting name of file
% d: dimension of vectors
% n: number of vectors
%
% OUTPUT
% F: matrix with unit-norm columns
%
% For more information, see https://www.math.colostate.edu/~king/GameofSloanes.html
% and the paper "Game of Sloanes: Best known packings in complex projective space"
% by John Jasper, Emily J. King, and Dustin G. Mixon

fid = fopen(fileName,'r');
F=fscanf(fid,'%f');

F=reshape(F(1:length(F)/2)+1i*F(length(F)/2+1:end),d,n);

F=F*diag(1./sqrt(diag(F'*F)));

fclose(fid);

end

