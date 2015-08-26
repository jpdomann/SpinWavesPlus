function [ out ] = reconstructMatrix( in,matrixSize )
%RECONSTRUCTMATRIX Input a column vector, output a 2D matrix

out = zeros(matrixSize);
for i = 1:numel(in)
    out(i) = in(i);
end

end

