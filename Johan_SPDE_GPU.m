function Johan_SPDE_GPU(sigma)
%%
% Path1 =  ;

%% Noise parameters 
% sigma = 0.1 ; % noise level
alpha = 0.1 ; % noise decay rate
kappa = 1 ; % related to time scale and noise

%% Spatial and temporal setting
N = 4096 ; dx = 0.5 ; dy = 0.5 ; Lengthx = N*dx ; Lengthy = N*dy ;
Grid = [N,N] ; Length = [Lengthx,Lengthy] ;
T = 1e2 + 0 ; dtref = 0.0005 ; dt = kappa*dtref;
Record = 1:2e5:T ; % Save data
Count = 1 ;

%% Model parameters
DifP = 0.1 ; DifW = 0.1 ; DifO = 100 ; % AdvO_X = 0 ; % 100*(rand(N)*2 - 1) ; AdvO_Y = 0 ; % 100*(rand(N)*2 - 1) ;
R = 0.7 ; alpha0 = 0.2 ; W0 = 0.2 ; rw = 0.2 ; cc = 10 ; 
gmax = 0.05 ; dd = 0.25 ; k1 = 5 ; k2 = 5 ;

%% Fourier setting
lambdax = 2*pi*[0:N/2,-N/2+1:-1]'/Lengthx ;
lambday = 2*pi*[0:N/2,-N/2+1:-1]'/Lengthy ;
[lambdaxx,lambdayy] = meshgrid(lambday,lambdax) ;
A = (lambdaxx.^2+lambdayy.^2) ;
MM_P = DifP*A ; EE_P = gpuArray(1./(1+dt*MM_P)) ;
MM_O = DifO*A ; EE_O = gpuArray(1./(1+dt*MM_O)) ;
MM_W = DifW*A ; EE_W = gpuArray(1./(1+dt*MM_W)) ;
bj = gpuArray(get_twod_bj(dtref,Grid,Length,alpha)) ; % get noise coeffs

%% Initial condition
P_pre = rand(N) ;
P = zeros(N) ; % Plant
P(P_pre < 0.05) = 100 ;
P = gpuArray(P) ;
O = gpuArray(zeros(N) + R/(alpha0*W0)) ; % Surface water
W = gpuArray(zeros(N) + R/rw/4) ; % Soil water
Ph = fft2(P) ; Oh = fft2(O) ; Wh = fft2(W) ;
tic
%% Loop
for ii = 1:T/kappa % time loop
    fPh = fft2(cc*gmax*W./(W + k1).*P - dd*P) ;
    fOh = fft2(R - alpha0*(P + k2*W0)./(P + k2).*O) ;
    fWh = fft2(alpha0*(P + k2*W0)./(P + k2).*O - gmax*W./(W + k1).*P - rw*W) ;
    [PdW,~] = get_twod_dW_GPU(bj,kappa) ;
    [OdW,~] = get_twod_dW_GPU(bj,kappa) ;
    [WdW,~] = get_twod_dW_GPU(bj,kappa) ;
    gPdWh = fft2(sigma.*PdW) ;
    gOdWh = fft2(sigma.*OdW) ;
    gWdWh = fft2(sigma.*WdW) ;
    Ph = EE_P.*(Ph+dt*fPh+gPdWh) ;
    Oh = EE_O.*(Oh+dt*fOh+gOdWh) ;
    Wh = EE_W.*(Wh+dt*fWh+gWdWh) ;
    P = real(ifft2(Ph)) ; O = real(ifft2(Oh)) ; W = real(ifft2(Wh)) ; 
%     if ismember(ii,Record) == 1
%         disp(['Time steps' num2str(ii)]) ;
%         DataW = gather(W) ; DataP = gather(P) ; DataO = gather(O) ;
%         Filename = ['Arid_R' num2str(R) '_sigma' num2str(sigma) '_' num2str(Count) '.mat'] ;
%         save([Path1 Filename],'DataP','DataW','DataO') ;
%         Count = Count + 1 ;
%     end
end
toc