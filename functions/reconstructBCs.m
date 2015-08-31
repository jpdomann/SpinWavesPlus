function [ final ] = reconstructBCs(initial,new,Indices,num, varargin )
%RECONSTRUCTBCS Takes the initially reconstructed data matrix, and
%reorients it to include the specified boundary conditions
%   Detailed explanation goes here


% Determine what row needs to be changed
nums = find(Indices == 1);
row = nums(num);


% Break initial data into two section (before and after desired row)
before = initial(1:row-1,:);
after = initial(row:end,:);

%Reconstruct
final = [before; new; after];

%Delete Last Row
if ~isempty(varargin)
    deleteLast = varargin{1};    
    if deleteLast
        final(end,:) = [];
    end
end

end

