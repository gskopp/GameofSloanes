File handling

convertToMat.m
Convert a directory filled with txt files of the required format to Matlab mat files named dxn_init.mat containing a dxn matrix with columns the vectors with variable name F.

convertToText.m
Convert a directory filled with Matlab mat files named dxn_init.mat containing a dxn matrix with columns the vectors with variable name F to a txt file of the required format. 

saveF.m
Given a dxn matrix with unit-norm columns, this file saves a txt file of the required format for submission as well as a mat file named dxn_init.mat containing the matrix with variable name F. 

saveText.m
Given a dxn matrix with unit-norm columns, this file saves a txt file of the required format for submission.

Implementations of algorithms from papers

GrassPack.m
Modification of the code from http://www.ece.mcmaster.ca/~davidson/pubs/Flexible_codebook_design.html which is the implementation of "Flexible codebook design for limited feedback systems via sequential smooth optimization on the Grassmannian manifold" by Ahmed Medra and Timothy N. Davidson

troppAltProj.m
Implementation of the algorithm in "Designing Structured Tight Frames Via an Alternating Projection Method" by Joel A. Tropp, Inderjit S. Dhillon, Robert W. Heath, Jr., and Thomas Strohmer and Tropp's dissertation "Topics in Sparse Approximation."