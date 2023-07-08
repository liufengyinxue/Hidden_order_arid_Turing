function Arid_Rietkerk_GPU(R)
%%
% This function simulate the model proposed by Rietkerk et al. (2002).
% Rietkerk, M., M. C. Boerlijst, F. van Langevelde, R. HilleRisLambers, J. van de Koppel, 
% L. Kumar, H. H. T. Prins, and A. M. de Roos. 2002. Self‐Organization of Vegetation in Arid Ecosystems. 
% The American Naturalist 160:524–530.

%%
GPU = gpuDevice(1);
reset(GPU) 

%% Spatial and temporal settings
N = 4096 ;
T = 6e6+1 ;
dx = 0.5 ;
dy = 0.5 ;
dt = 0.0005 ;

%% Parameters
% R = 0.75 ;
DifP = 0.1 ; DifW = 0.1 ; DifO = 100 ;
alpha = 0.2 ; W0 = 0.2 ; rw = 0.2 ;
cc = 10 ; gmax = 0.05 ; dd = 0.25 ;
k1 = 5 ; k2 = 5 ;

%% Kernel
kernel = gpuArray([0,1/(dx*dx),0;1/(dy*dy),-2/(dx*dx)-2/(dy*dy),1/(dy*dy);0,1/(dx*dx),0]) ;

%% Initial condition
P_pre = rand(N) ;
P = zeros(N) ; % Plant
P(P_pre < 0.05) = 100 ;
P = gpuArray(P) ;
% P(P_pre >= 0.05) = 10 ;
O = gpuArray(zeros(N) + R/(alpha*W0)) ; % Surface water
W = gpuArray(zeros(N) + R/rw/4) ; % Soil water

%% save data
Record = 1:6e4:T ;
Count = 1 ;

%% Loop
for i = 1:T
    d2Odxy2 = DifO*imfilter(O,kernel,'circular') ;
    drO = R - alpha*(P + k2*W0)./(P + k2).*O ;
    d2Wdxy2 = DifW*imfilter(W,kernel,'circular') ;
    drW = alpha*(P + k2*W0)./(P + k2).*O - gmax*W./(W + k1).*P - rw*W ;
    d2Pdxy2 = DifP*imfilter(P,kernel,'circular') ;
    drP = cc*gmax*W./(W + k1).*P - dd*P ;
    O = O + (drO + d2Odxy2)*dt ;
	W = W + (drW + d2Wdxy2)*dt ;
	P = P + (drP + d2Pdxy2)*dt ;
%     if ismember(i,Record) == 1
%         disp(['Time steps' num2str(i)]) ;
%         DataP = gather(P) ; DataO = gather(O) ; DataW = gather(W) ;
%         Filename = ['Arid_R' num2str(R) '_' num2str(Count) '.mat'] ;
%         save([Path Filename],'DataP','DataO','DataW') ;
%         Count = Count + 1 ;
%     end
end