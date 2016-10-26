function [SBMLNotes] = convert2SBMLNotes(SBmodelNotes)
% convert2SBMLNotes
% converts the notes string given by SBTOOLBOX that they can be used
% through the SBML Toolbox (conversion to XHTML)
%
% USAGE:
% ======
% SBMLNotes = convert2SBMLNotes(SBmodelNotes)
%
% SBmodelNotesString: a notes string created manually at the commmand
%                     window or in one of the GUIs
% SBMLNotes: a string in XHTML format that can be converted by SBML Toolbox

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% programmed by Gunnar Drews, Student at Department of Computerscience
% University of Rostock, gunnar.drews@uni-rostock.de
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

SBMLNotes = '';
newline = sprintf('\n');
notesStart = char([double('<html xmlns="http://www.w3.org/1999/xhtml">'), double(newline), double('<body>'), double(newline)]);
notesEnd = char([double('</body>'), double(newline), double('</html>')]);

% test wether input string is empty
if (length(SBmodelNotes)~=0),
    % test wether notes do already contain xhtml tags
   if ((length(strfind(SBmodelNotes, '<html'))~=0) || (length(strfind(SBmodelNotes, '<body'))~=0) || (length(strfind(SBmodelNotes, '<p>'))~=0)), 
       warn = sprintf('WARNING: "SBmodelNotes" already contain XHTML tags, conversion aborted!\nNotes are passed through without changes.')
       SBMLNotes = SBmodelNotes;
   else
       % test wether notes field contains several lines
       linebreaks = strfind(SBmodelNotes, newline);
       if (length(linebreaks)==0),
           % build single line SBML notes
           SBMLNotes = char([double(notesStart), double(SBmodelNotes), double(notesEnd)]);
       else
           % build multiple line SBML notes
           SBMLNotes = notesStart;
           startPos = 1;
           for index = 1 : length(linebreaks),
               endPos = linebreaks(index)-1;
               SBMLNotes = char([double(SBMLNotes), double('<p>'), double(SBmodelNotes(startPos:endPos)), double('</p>'), double(newline)]);
               startPos = linebreaks(index)+1;
           end
           SBMLNotes = char([double(SBMLNotes), double('<p>'), double(SBmodelNotes(startPos:length(SBmodelNotes))), double('</p>'), double(newline), double(notesEnd)]);
       end
   end
end
return