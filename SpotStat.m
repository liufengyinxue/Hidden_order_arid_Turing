function [SpotNum, SpotNum_density] = SpotStat(DataP)
%%% This function compute the number of spot vegtation patches in spatial
%%% patterns.

    Circle_Thesh = 0.9 ;
    A = DataP ;
    B = imbinarize(A) ;
    [C,Num] = bwlabel(B,8) ;
    Obj = regionprops(C,'Circularity') ;
    Circularity = cat(1,Obj.Circularity) ; 
    [Index1,~] = find(Circularity < Circle_Thesh) ;
%     Ctemp = C ;
%     Ctemp(ismember(C,Index1)) = 0 ;
%     D = Ctemp > 0 ;
    SpotNum = Num - length(Index1) ;
    SpotNum_density = SpotNum/mean(A,'all') ;
end