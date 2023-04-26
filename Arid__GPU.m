function Arid_Johan_GPU(R)
GPU = gpuDevice(1);
reset(GPU) 
Path = 'F:\ZhenpengGe\Arid_data\' ;
N = 4096 ;
T = 6e7+1 ;

DifP = 0.1 ;
DifW = 0.1 ;
DifO = 100 ;
AdvO = 0 ;
% R = 0.75 ;
alpha = 0.2 ;
W0 = 0.2 ;
rw = 0.2 ;
cc = 10 ;
gmax = 0.05 ;
dd = 0.25 ;
k1 = 5 ;
k2 = 5 ;

dx = 0.5 ;
dy = 0.5 ;
dt = 0.0005 ;
kernel = gpuArray([0,1/(dx*dx),0;1/(dy*dy),-2/(dx*dx)-2/(dy*dy),1/(dy*dy);0,1/(dx*dx),0]) ;
kernelx = gpuArray([0,0,0;-1/(2*dx),0,1/(2*dx);0,0,0]) ;
kernely = gpuArray([0,-1/(2*dy),0;0,0,0;0,1/(2*dy),0]) ;

P_pre = rand(N) ;
P = zeros(N) ; % Plant
P(P_pre < 0.05) = 100 ;
P = gpuArray(P) ;
% P(P_pre >= 0.05) = 10 ;
O = gpuArray(zeros(N) + R/(alpha*W0)) ; % Surface water
W = gpuArray(zeros(N) + R/rw/4) ; % Soil water

% Record1 = 1:100:1e5 ;
% Record2 = ceil(logspace(5.1,6,100)) ;
Record = 1:6e4:T ;
Count = 1 ;

F1 = figure('position',[200 20 600 600]) ;
Videoname = ['Johan_arid_' num2str(N) '_T' num2str(T)] ;
V1 = VideoWriter(Videoname,'MPEG-4') ; % 'MPEG-4'
V1.FrameRate = 30 ;
V1.Quality = 100 ;
open(V1) ;
for i = 1:T
    if rem(i - 1,1e4) == 0
        figure(F1) ;
        imagesc(gather(P)) ;
        axis equal
        xlim([0 N]) ;
        ylim([0 N]) ;
        colorbar ; 
%         caxis([0,1]) ;
        title(['Plant density, T=' num2str(i)]) ;
        Frame = getframe(gcf) ;
        writeVideo(V1,Frame) ;
    end
    d2Odxy2 = DifO*imfilter(O,kernel,'circular') ; % - AdvO*imfilter(O,kernely,'circular') ;
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
%         DataP = gather(P) ;
%         Filename = ['Arid_R' num2str(10*R) '_' num2str(Count) '.mat'] ;
%         save([Path Filename],'DataP') ;
%         Count = Count + 1 ;
%     end
%     if i == T
%         DataP = gather(P) ;
%         DataO = gather(O) ;
%         DataW = gather(W) ;
%         Filename_final = ['Arid_R' num2str(10*R) '_T6e6plus1_final.mat'] ;
%         save([Path Filename_final],'DataP','DataO','DataW') ;
%     end
end
close(V1) ;