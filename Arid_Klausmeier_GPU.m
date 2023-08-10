function Arid_Klausmeier_GPU(a)
%%
% This function simulate the model proposed by Klausmeier (1999).
% Klausmeier, C. A. 1999. Regular and Irregular Patterns in Semiarid Vegetation. 
% Science 284:1826–1828.

%%
GPU = gpuDevice(1);
reset(GPU) 

%% Spatial and temporal settings
N = 4096 ;
T = 2e7 + 1 ;
dx = 0.5 ;
dy = 0.5 ;
dt = 0.0001 ;

%% Parameters
% a = 0.6 ;
m = 0.45 ;
Dw = 500 ;

%% Kernel
kernel = gpuArray([0,1/(dx*dx),0;1/(dy*dy),-2/(dx*dx)-2/(dy*dy),1/(dy*dy);0,1/(dx*dx),0]) ;

%% Initial condition
w = 2*ones(N) ; 
nn = rand(N) ;
n = zeros(N) ;
n(nn < 0.1) = 10 ;
w = gpuArray(w) ; n = gpuArray(n) ;

%% save data
% Record = 1:2e5:T ;
% Count = 1 ;

%% Loop
for i = 1:T
    wi = a - w - w.*n.^2 + Dw*imfilter(w,kernel,'circular') ; % V*imfilter(w,kernely,'circular'), 'symmetric' 'circular'
    ni = w.*n.^2 - m*n + imfilter(n,kernel,'circular') ;
    w = w + dt*wi ;
    n = n + dt*ni ;
%     if ismember(i,Record) == 1
%         disp(['Time steps' num2str(i)]) ;
%         DataW = gather(w) ; DataP = gather(n) ;
%         Filename = ['Arid_a' num2str(a) '_' num2str(Count) '.mat'] ;
%         save([Path Filename],'DataP','DataW') ;
%         Count = Count + 1 ;
%     end
end