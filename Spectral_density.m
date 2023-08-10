function [SK] = Spectral_density(Data) 
%%%
% This function computes circularly averaged spectral density
% for 2D matrix data.
% Input:     
% Data:    a 2D matrix, like an image
% Output:  
% SK:      a  matrix, the first column is wave number,
%             the second column is circularly averaged spectral density.
% % % 

% mat2gray
img = floor(255*mat2gray(Data));
img = 255 - uint8(img);
mean_gray = mean(img(:));
img = img - mean_gray;
% perform 2D fft, compute psd
[N,M] = size(img);
% Compute power spectrum % with gpu
img = single(img);
% img = gpuArray(img);
imgf = fftshift(fft2(img));
% imgf = gather(imgf);
% img = gather(img);
imgfp = (abs(imgf)/(N*M)).^2;   
%%% Adjust PSD size
dimDiff = abs(N-M);
dimMax = max(N,M);
if N > M                                                                    % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
    end
elseif N < M                                                                % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
end
% halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)
%%% Compute radially average power spectrum
[X,Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[theta,rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
rho = round(rho)+1;
rho = rho(:);
imgfp= imgfp(:);
Pf = zeros(max(rho),1);
Counter = Pf;
for index = 1:length(rho)
    Pf(rho(index)) = Pf(rho(index)) + imgfp(index); 
    Counter(rho(index)) = Counter(rho(index)) + 1; 
end
Pf = Pf./Counter; % circular average
Pf(isnan(Pf)) = 0;
res = 0.5 ;
Fs = 1/2*(1/res);
L = length(Pf)-1;
f = Fs*(0:L)'/L;
SK=[f Pf];
SK(1,:)=[];
SK(SK(:,1)>0.3*(1/res),:)=[];
SK(SK(:,1)<4/size(img,1)*(1/res),:)=[];
end