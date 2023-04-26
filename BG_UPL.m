load( 'Arid_R0.9_20.mat') ; % load simulation data
[M1,N1] = size(DataP) ;
dx = 0.5 ; % grid size, unit: m
Threshold = 1 ;
BW = imbinarize(DataP,Threshold) ; % Divide DataP into vegetation (1) and bare ground (0)
[Index1, Index2] = find(BW == 0) ; % find the index of bare ground
% below: calculate 8-neighborhood gradient
% North = [1;0;0] ; South = [0;1] ; East = [0,1] ; West = [1,0,0] ;
N = [0,1,0;0,0,0;0,0,0] ; S = [0,0,0;0,0,0;0,1,0] ; 
E = [0,0,0;0,0,1;0,0,0] ; W = [0,0,0;1,0,0;0,0,0] ;
NE = [0,0,1;0,0,0;0,0,0] ; NW = [1,0,0;0,0,0;0,0,0] ; 
SE = [0,0,0;0,0,0;0,0,1] ; SW = [0,0,0;0,0,0;1,0,0] ; 
O_N = imfilter(DataO,N,'circular') ; O_S = imfilter(DataO,S,'circular') ;
O_E = imfilter(DataO,E,'circular') ; O_W = imfilter(DataO,W,'circular') ;
O_NE = imfilter(DataO,NE,'circular') ; O_NW = imfilter(DataO,NW,'circular') ;
O_SE = imfilter(DataO,SE,'circular') ; O_SW = imfilter(DataO,SW,'circular') ;
Slope_N = (DataO - O_N)/dx ; Slope_S = (DataO - O_S)/dx ;
Slope_E = (DataO - O_E)/dx ; Slope_W = (DataO - O_W)/dx ;
Slope_NE = (DataO - O_NE)/(sqrt(2)*dx) ; Slope_NW = (DataO - O_NW)/(sqrt(2)*dx) ;
Slope_SE = (DataO - O_SE)/(sqrt(2)*dx) ; Slope_SW = (DataO - O_SW)/(sqrt(2)*dx) ;
% 
State_N = imfilter(BW,N,'circular') ; State_S = imfilter(BW,S,'circular') ;
State_E = imfilter(BW,E,'circular') ; State_W = imfilter(BW,W,'circular') ;
State_NE = imfilter(BW,NE,'circular') ; State_NW = imfilter(BW,NW,'circular') ;
State_SE = imfilter(BW,SE,'circular') ; State_SW = imfilter(BW,SW,'circular') ;
% below: we don't calcuate UPL of bourdary points
Abandon_length = 100 ;
Abandon_row_index = [1:Abandon_length,M1+1-Abandon_length:M1] ;
Abandon_column_index = [1:Abandon_length,M1+1-Abandon_length:M1] ;
%
UPL = zeros(length(Index1),1) ;
LabelNum = 1e6 ; % A label represents the lowest point and water cannot reach the vegetation patch
for ii = 1:length(Index1)
    if ismember(Index1(ii),Abandon_row_index) == 1 && ismember(Index2(ii),Abandon_column_index) == 1
        UPL(ii) = nan ;
    else
%         Sign = zeros(M1,N1) ; % it's not necessary for calculating UPL
        IndX = Index1(ii) ; IndY = Index2(ii) ;
%         Sign(IndX,IndY) = 1 ;
        n1 = 0 ; n2 = 0 ;
        % identify whether the neighbor is a vegetation patch 
        while (State_N(IndX,IndY) == 0) && (State_S(IndX,IndY) == 0) && (State_E(IndX,IndY) == 0)  && (State_W(IndX,IndY) == 0) && (State_NE(IndX,IndY) == 0) && (State_NW(IndX,IndY) == 0) && (State_SE(IndX,IndY) == 0) && (State_SW(IndX,IndY) == 0)
            Slope = [Slope_N(IndX,IndY),Slope_S(IndX,IndY)...
                Slope_E(IndX,IndY),Slope_W(IndX,IndY)...
                Slope_NE(IndX,IndY),Slope_NW(IndX,IndY)...
                Slope_SE(IndX,IndY),Slope_SW(IndX,IndY)] ;
            if sum(double(Slope < 0)) == 8 % lowest point
                n1 = LabelNum/dx/2 ; n2 = LabelNum/(sqrt(2)*dx)/2 ; % just a label to present that surface water can not reach vegetation patch
                break
            else
                [~,Position] = max(Slope) ;
                if Position == 1
                    IndX = IndX - 1 ; n1 = n1 + 1 ; % North
                elseif Position == 2
                    IndX = IndX + 1 ; n1 = n1 + 1 ; % South
                elseif Position == 3
                    IndY = IndY + 1 ; n1 = n1 + 1 ; % East
                elseif Position == 4
                    IndY = IndY - 1 ; n1 = n1 + 1 ; % West
                elseif Position == 5
                    IndX = IndX - 1 ; IndY = IndY + 1 ; n2 = n2 + 1 ; % NorthEast
                elseif Position == 6
                    IndX = IndX - 1 ; IndY = IndY - 1 ; n2 = n2 + 1 ; % NorthWest
                elseif Position == 7
                    IndX = IndX + 1 ; IndY = IndY + 1 ; n2 = n2 + 1 ; % SouthEast
                elseif Position == 8
                    IndX = IndX + 1 ; IndY = IndY - 1 ; n2 = n2 + 1 ; % SouthWest
                end
            end
%             Sign(IndX,IndY) = 1 ;
        end
        UPL(ii) = n1*dx + n2*sqrt(2)*dx ;
    end
end
% Calculate spatial mean UPL
UPL(UPL == LabelNum) = nan ;
mUPL = mean(UPL,'omitnan') ;