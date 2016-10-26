% buildMEXfilesWINDOWS: Builds the lipsol mex files on Windows systems
% In order to build the MEX files on Windows you need to install the MinGW
% compiler.
%
% Please note that usually on Windows the MEX files do not have to be rebuild.

% build the MEX files using the mexfSB function included in the SBTOOLBOX
cd srcwindows
mexfSB('ordmmd.f ordmmdg.f');
mexfSB('symfct.f  symfctg.f');
mexfSB('inpnv.f   inpnvg.f');
mexfSB('blkfct.f  blkfctg.f');
mexfSB('blkslv.f  blkslvg.f');

% clean up
!move *.mexw32 ../mex/
!delete *.obj
cd ..
clc
disp('lipsol MEX files ready!');

