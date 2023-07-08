function [VMR] = Raster_DF_ThrowBox(Data,Scale,BoxNum,Pbc)
% This function computes density fluctuation
% of raster data based on the method of random boxes (square box)¡£
% Input:
% Data:       a matrix, both of binary data (consisting of '1' and '0') and
%                 scalar field are fine.
% Scale:      the sizes of the boxes, must be a column vector
% Pbc:         Indicates whether periodic boundary conditions are adopted or not, 
%                 generally, it is necessary to adopt for simulation data with periodic
%                 boundary conditions, turn on: 1 (periodic boundary condition), 
%                 turn off: 0 (open boundary condition).
% BoxNum: the number of times of throwing boxes, a number
% Output:
% VMR:       the density fluctuation

v = zeros(length(Scale),BoxNum) ;
switch Pbc
    case 0
        [M1,N1] = size(Data) ;
        for i = 1:length(Scale)
            for j = 1:BoxNum
                Index1 = randi([1,M1 - Scale(i) + 1]) ; %the row index of top-left corner
                Index2 = randi([1,N1 - Scale(i) + 1]) ; %the column index of top-left corner
                Indexrow = Index1:1:(Index1+Scale(i)-1) ;
                Indexcol = Index2:1:(Index2+Scale(i)-1) ;
                tempdata = Data(Indexrow,Indexcol) ;
                v(i,j) = sum(tempdata(:)) ;
            end
        end
    case 1
        for i = 1:length(Scale)
            padsize = Scale(i) - 1 ;
            ExpandedData = padarray(Data,[padsize padsize],'circular') ;
            [M2,N2] = size(ExpandedData) ;
            for j = 1:BoxNum
                Index1 = randi([1,M2 - Scale(i) + 1]) ; %the row index of top-left corner
                Index2 = randi([1,N2 - Scale(i) + 1]) ; %the column index of top-left corner
                Indexrow = Index1:1:(Index1+Scale(i)-1) ;
                Indexcol = Index2:1:(Index2+Scale(i)-1) ;
                tempdata = ExpandedData(Indexrow,Indexcol) ;
                v(i,j) = sum(tempdata(:)) ;
            end
        end
    otherwise
end
Square_scales = repmat(Scale,1,BoxNum) ;
VMR = var(v./Square_scales.^2,0,2) ;