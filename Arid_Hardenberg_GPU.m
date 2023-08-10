function Arid_Hardenberg_GPU(p)
%%
% This function simulate the model proposed by von Hardenberg et al. (2002).
% von Hardenberg, J., E. Meron, M. Shachak, and Y. Zarmi. 2001. 
% Diversity of Vegetation Patterns and Desertification. 
% Physical Review Letters 87:198101.

%% Spatial and temporal settings
N = 4096 ;
T = 1e7 + 1 ;
dx = 0.5 ; % 1
dy = 0.5 ; % 1
dt = 0.0005 ; % 0.002

%% Parameters
gamma = 1.6 ;
sigma = 1.6 ;
mu = 0.2 ;
pho = 1.5 ;
delta = 100 ;
% p = 0.45 ; % 0-0.6
beta = 3 ;

%% Kernel
kernel = gpuArray([0,1/(dx*dx),0;1/(dy*dy),-2/(dx*dx)-2/(dy*dy),1/(dy*dy);0,1/(dx*dx),0]) ;

%% Initial condition
% % a>2*m
% n = (a + sqrt(a^2 - 4*m^2))/(2*m) + 0.01*(2*rand(N) - 1) ;
% w = 2*m^2/(a + sqrt(a^2 - 4*m^2)) + zeros(N) ;
w = 1 + zeros(N) ; % p
nn = rand(N) ;
n = zeros(N) ;
n(nn < 0.1) = 1 ;
w = gpuArray(w) ; n = gpuArray(n) ;

%% save data
Record = 1:2e5:T ;
Count = 1 ;

%% Loop
for i = 1:T
    wi = p - (1 - pho*n).*w - w.^2.*n + delta*imfilter(w - beta*n,kernel,'circular') ; % 'symmetric' 'circular'
    ni = (gamma*w)./(1 + sigma*w).*n - n.^2 - mu*n + imfilter(n,kernel,'circular') ;
    w = w + dt*wi ; n = n + dt*ni ;
%     if ismember(i,Record) == 1
%         disp(['Time steps' num2str(i)]) ;
%         DataW = gather(w) ; DataP = gather(n) ;
%         Filename = ['Arid_p' num2str(p) '_' num2str(Count) '.mat'] ;
%         save([Path Filename],'DataP','DataW') ;
%         Count = Count + 1 ;
%     end
end