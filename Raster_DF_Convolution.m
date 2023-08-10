function [VMR] = Raster_DF_Convolution(Data,Scale,BoxNum,Pbc,GPU)
% This function computes density flucutation
% based on the method of convolution (square box).
% This function can be executed on GPU.
% Input:
% Data:       a matrix
% Scale:      the sizes of the boxes, a vector
% Pbc:         Indicates whether periodic boundary conditions are adopted or not, 
%                 generally, it is necessary to adopt for simulation data with periodic
%                 boundary conditions, turn on: 1, turn off: 0
% GPU:       Indicates whether it is executed on the GPU, turn on: 1, turn off: 0
% BoxNum: the number of Boxes
% Output:
% VMR:       the density fluctuation

[L1,L2] = size(Data) ;
if GPU == 1
    DataTemp = gpuArray(Data) ;
elseif GPU == 0 
    DataTemp = Data ;
end

VMR = zeros(length(Scale),1); 
for i = 1:length(Scale)
    kernel = ones(Scale(i)) ; 
    kernel = single(kernel) ; 
    kersum = sum(kernel(:)) ;
    if GPU == 1
        kerneltemp = gpuArray(kernel) ; 
    elseif GPU == 0
        kerneltemp = kernel ;
    end
    if Pbc == 0
        DesityTemp = conv2(DataTemp,kerneltemp,'valid') ;
        if GPU == 1
            DesityTemp_C = gather(DesityTemp) ;
            DesityTemp_C_Flatten = DesityTemp_C(:) ;
            SelectBox = randi((L1 - Scale(i) + 1)*(L2 - Scale(i) + 1),[BoxNum,1]) ;
            VMR(i,1) = var(DesityTemp_C_Flatten(SelectBox)/kersum) ;
        elseif GPU == 0
            DesityTemp_Flatten = DesityTemp(:) ;
            SelectBox = randi((L1 - Scale(i) + 1)*(L2 - Scale(i) + 1),[BoxNum,1]) ;
            VMR(i,1) = var(DesityTemp_Flatten(SelectBox)/kersum) ;
        end
    elseif Pbc == 1
        padsize = Scale(i) - 1 ;
        DT = padarray(DataTemp,[padsize padsize],'circular') ; 
        DesityTemp = conv2(DT,kerneltemp,'valid') ;
        if GPU == 1
            DesityTemp_C = gather(DesityTemp) ;
            DesityTemp_C_Flatten = DesityTemp_C(:) ;
            SelectBox = randi(L1*L2,[BoxNum,1]) ;
            VMR(i,1) = var(DesityTemp_C_Flatten(SelectBox)/kersum) ;
        elseif GPU == 0
            DesityTemp_Flatten = DesityTemp(:) ;
            SelectBox = randi(L1*L2,[BoxNum,1]) ;
            VMR(i,1) = var(DesityTemp_Flatten(SelectBox)/kersum) ;
        end
    end
end
end