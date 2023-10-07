function [data, sigma] = length_norm(data)
% length-normalizes the data vectors assuming one observation per column.
%
% Inputs:
%   - data        : input data matrix, one observation per column
% Outputs:
%   - data        : output length-normalized data matrix     
%
% Omid Sadjadi <s.omid.sadjadi@gmail.com>
% Microsoft Research, Conversational Systems Research Center

sigma = sqrt(sum(data.^2));
data  = bsxfun(@rdivide, data, sigma);
