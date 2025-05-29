function saveText( F, init )
%saveText 
% Save packing in text format a la Sloane, as a single column,
% first of the real parts then the imaginary.
%
% INPUTS:
% F : a matrix with unit-norm columns
% init : initials of the creator of the packing

[d n]=size(F);

F=F(:);
F=vertcat(real(F),imag(F));

fid = fopen(strcat(num2str(d),'x',num2str(n),'_',init,'.txt'),'w');
for ii = 1:2*d*n
    fprintf(fid,'%.15f\n',F(ii));
end
fclose(fid);


end

