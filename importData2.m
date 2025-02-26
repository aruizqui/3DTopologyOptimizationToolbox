function [elemConn, nodeCoord] = importData2(filename)

% Import data from .fem file
% Creates text file coordinates.txt and connectivity.txt using python
% and imports the data from these files
% v1 Feb 2021 - V. Cholvi: Inital version Toolbox
% v2 Dec 2024 - A. Ruiz: update root and clean code

python_script = 'TopologyOptimizationToolbox/transformData_v2.py';
command = sprintf('python "%s" "%s"',python_script, filename);
system(command);

nodeCoord = importdata('coordinates.txt');

elemConn = importdata('conectivity.txt');

checkInvalid = sum(isnan([nodeCoord(:);elemConn(:)])) ;
if (checkInvalid > 0) 
    error('__Invalid Data in Imported Arrays')
end

end
