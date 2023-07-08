function Arid_Getzin_GPU(p)
%%
% This function simulate the model proposed by Getzin et al. (2016).
% Getzin, S., H. Yizhaq, B. Bell, T. E. Erickson, A. C. Postle, I. Katra, O. Tzuk, Y. R. Zelnik, K. Wiegand, T. Wiegand, and E. Meron. 2016. 
% Discovery of fairy circles in Australia supports self-organization theory. Proceedings of the National Academy of Sciences 113:3551–3556.

%%
GPU = gpuDevice(1);
reset(GPU) 

%% Spatial and temporal settings
N = 4096 ;
T = 2e7 + 1 ;
dx = 1 ; 
dy = 1 ; 
dt = 0.00004 ;

%% Parameters
% P = 280 ;
K = 0.666 ; Q = 1.2 ; M = 2 ; A = 120 ; NW = 1.5 ; NH = 4.5 ; 
E = 1.5 ; Lambda = 0.03 ; Gamma = 14 ; DB = 0.1 ; DW = 2.5 ; 
DH = 4 ;  RW = 0.3 ; RH = 0.8 ; f = 0.01 ;

% p = P*Lambda/(M^2) ; 
q = Q/K ; nuw = NW/M ; nuh = NH/M ; alpha = A/M ; eta = E*K ;
gamma = Gamma*K/M ; deltaw = DW/DB ; deltah = DH*M/(DB*Lambda) ;
Rw = RW ; Rh = RH ;

%% Kernel
kernel = gpuArray([0,1/(dx*dx),0;1/(dy*dy),-2/(dx*dx)-2/(dy*dy),1/(dy*dy);0,1/(dx*dx),0]) ;

%% Initial condition
w = alpha*f*p/(nuw*(alpha*f + nuh)) + zeros(N) ; 
h = p/(alpha*f + nuh) + zeros(N) ;
bb = rand(N) ;
b = zeros(N) ;
b(bb < 0.1) = 1 ;
w = gpuArray(w) ; h = gpuArray(h) ; b = gpuArray(b) ;

%% save data
Record = 1:2e5:T ;
Count = 1 ;

%% Loop
for i = 1:T
    wi = alpha*(b + q*f)./(b + q).*h - nuw./(1 + Rw*b).*w - gamma*b.*(1 + eta*b).^2.*w + deltaw*imfilter(w,kernel,'circular') ;
    hi = p - alpha*(b + q*f)./(b + q).*h - nuh./(1 + Rh*b).*h + deltah*imfilter(h.^2,kernel,'circular') ;
    bi = w.*(1 + eta*b).^2.*b.*(1 - b) - b + imfilter(b,kernel,'circular') ;
    w = w + dt*wi ; h = h + dt*hi ; b = b + dt*bi ;
%     if ismember(i,Record) == 1
%         disp(['Time steps' num2str(i)]) ;
%         DataW = gather(w) ; DataO = gather(h) ; DataP = gather(b) ;
%         Filename = ['Arid_p' num2str(p) '_' num2str(Count) '.mat'] ;
%         save([Path Filename],'DataP','DataW','DataO') ;
%         Count = Count + 1 ;
%     end
end