%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/yanxingsu/Documents/USC/ASTE 585 spacecraft control/horizons_results.txt
%
% Auto-generated by MATLAB on 21-Jul-2024 15:12:00

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 18);

% Specify range and delimiter
opts.DataLines = [1020, 1824];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "updated", "Mar04", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"];
opts.SelectedVariableNames = ["updated", "Mar04"];
opts.VariableTypes = ["string", "string", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "updated", "TrimNonNumeric", true);
opts = setvaropts(opts, "updated", "ThousandsSeparator", ",");

% Import the data
horizonsresults1 = readtable("/Users/yanxingsu/Documents/USC/ASTE 585 spacecraft control/horizons_results.txt", opts);


%% Clear temporary variables
clear opts