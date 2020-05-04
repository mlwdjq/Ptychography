
% Add this root directory to path
[cDirThis, cName, cExt] = fileparts(mfilename('fullpath'));
addpath(genpath(cDirThis));
addpath('../mpm');

% mic library
mpm addpath


% Build PIE UI
pie = PIE.ui.PIE_Analyze;


pie.build(-100, -50);

% If ryan's computer, then put on leftmost monitor for now
if strcmp(char(java.lang.System.getProperty('user.name')), 'T470P')
    drawnow
%    lsi.hFigure.Position = [-3966         500        1600         1000];
%    lsi.hFigure.Position = [-3966         500        1600         1000];
end