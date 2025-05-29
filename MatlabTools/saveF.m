%saveF 
% Code to save a packing as both a Matlab and text file with the 
% desired filename format. The data is written to a text file in a single 
% column, first of the real components of each vector then the imaginary.
%
% INPUTS:
% F = matrix with unit-norm columns
% init = string of three-four characters, labeling the sbumission
%
% For more information, see https://www.math.colostate.edu/~king/GameofSloanes.html
% and the paper "Game of Sloanes: Best known packings in complex projective space"
% by John Jasper, Emily J. King, and Dustin G. Mixon

function saveF(F,init)

[d n]=size(F);

save(strcat(num2str(d),'x',num2str(n),'_',init,'.mat'),'F')

F=F(:);
F=vertcat(real(F),imag(F));
fid = fopen(strcat(num2str(d),'x',num2str(n),'_',init,'.txt'),'w');
for ii = 1:2*d*n
    fprintf(fid,'%.15f\n',F(ii));
end
fclose(fid);

