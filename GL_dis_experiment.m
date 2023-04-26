function GL_dis_experiment(R,Proportion,Time)
GPU = gpuDevice(1);
reset(GPU) 
Path1 = ['G:\Arid_Johan_data\Longterm_simulation\Arid_Johan_R' num2str(R) '_data_longterm\'] ;
Filename1 = ['Arid_R' num2str(R) '_' num2str(Time) '.mat'] ;
Full_Data = load([Path1 Filename1]) ;
DataP1 = Full_Data.DataP ; DataO = Full_Data.DataO ; DataW = Full_Data.DataW ; 
Initial_Biomass = mean(DataP1,'all') ;
[DataP2] = GL_dis(DataP1,Proportion) ;
Biomass = mean(DataP2,'all') ;
P = gpuArray(DataP2) ; O = gpuArray(DataO) ; W = gpuArray(DataW) ;

DifP = 0.1 ;
DifW = 0.1 ;
DifO = 100 ;
% AdvO = 0 ;
% R = 0.70 ;
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

Count = 0 ; 
while Biomass < Initial_Biomass
    Count = Count + 1 ;
    d2Odxy2 = DifO*imfilter(O,kernel,'circular') ; % - AdvO*imfilter(O,kernely,'circular') ;
    drO = R - alpha*(P + k2*W0)./(P + k2).*O ;
    d2Wdxy2 = DifW*imfilter(W,kernel,'circular') ;
    drW = alpha*(P + k2*W0)./(P + k2).*O - gmax*W./(W + k1).*P - rw*W ;
    d2Pdxy2 = DifP*imfilter(P,kernel,'circular') ;
    drP = cc*gmax*W./(W + k1).*P - dd*P ;
    O = O + (drO + d2Odxy2)*dt ;
	W = W + (drW + d2Wdxy2)*dt ;
	P = P + (drP + d2Pdxy2)*dt ;
    DataP3 = gather(P) ; Biomass = mean(DataP3,'all') ;
	if rem(Count,1e4) == 0
        disp(['R' num2str(R) ' Time:' num2str(Time) ' LossProportion' num2str(Proportion) ' Count:' num2str(Count) ' Biomass:' num2str(Biomass)]) ;
    end
end
disp(['Completed Time:' num2str(Time) ' Count:' num2str(Count)]) ;
DataP = gather(P) ; DataO = gather(O) ; DataW = gather(W) ;
Path2 = ['G:\Arid_Johan_data\Disturbance_experiment\R' num2str(R) '_GL_vary_Proportion\'] ;
Filename2 = ['Arid_R' num2str(R) '_' num2str(Time) '_LossProportion' num2str(Proportion) '_GL_Recovery.mat'] ;
save([Path2 Filename2],'DataP','DataO','DataW','DataP2','Count') ;

function [BBB] = GL_dis(AAA,Proportion)
BBB = AAA*(1 - Proportion) ;
end
end