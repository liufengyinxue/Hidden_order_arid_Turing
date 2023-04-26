clear all ; close all ; clc ;

Path1 = 'C:\Workspace\WritingPapers\No.27Turing_hyperuniform\Data\Arid_Johan_data\' ;
Filename1 = ['Arid_R9_' num2str(10) '.mat'] ;
load([Path1 Filename1]) ;
[M,N] = size(DataP) ;
dx = 0.5 ;
Threshold = mean(DataP(:)) ; 
Data = double(imbinarize(DataP,Threshold)) ;
imagesc(Data) ; colorbar ;
Points_number = 1 ; % 1e4
Nearest_distance = zeros(Points_number,1) ; % "0" means selecting plant points. 
Length = 100  ;
sq = ones(3) ;
for i1 = 1:Points_number
    % choose a random point
    Choose_point = [randi([Length + 1,M - Length]),randi([Length + 1,N - Length])] ; % randi([Length + 1,M - Length])
    if Data(Choose_point(1),Choose_point(2)) == 0
        % construct the room of interest
        Spatial_extent = [Choose_point(1) - Length, Choose_point(1) + Length, Choose_point(2) - Length, Choose_point(2) + Length] ;
        ROI = zeros(M,N) ;
        ROI(Spatial_extent(1):Spatial_extent(2),Spatial_extent(3):Spatial_extent(4)) = 1 ;
        Data_temp = Data.*ROI ;
        PostBIm = Data_temp > 0 ; 
        [PostL,PostNum] = bwlabel(PostBIm,8) ;
        Temp = zeros(PostNum,1) ;
        for i2 = 1:PostNum
            ClusterTemp = PostL == i2 ;
            EB = imerode(ClusterTemp,sq) ;
            Cluster_int = ClusterTemp & ~EB ;
            [x1,y1] = find(Cluster_int == 1) ;
            Cord1 = [x1,y1] ;
            Temp(i2,1) = min(pdist2(Choose_point,Cord1)) ;
        end
        Nearest_distance(i1,1) = min(Temp)*dx ;
    else
        continue ;
    end
end

%% check result
figure(1) ; imagesc(Data) ; impixelinfo ;
figure(2) ; imagesc(Data_temp) ; impixelinfo ;