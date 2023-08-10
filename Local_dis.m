function [BBB] = Local_dis(AAA,Proportion)
    Loss_Biomass = mean(AAA,'all')*Proportion ;
    AA = rescale(AAA) ;
    B = imbinarize(AA,'adaptive') ; % 'adaptive'
    [C,Patch_Num] = bwlabel(B,4) ;
    P_Remove_Patch = randperm(Patch_Num) ;
    Patch_Mean_Biomass = zeros(length(P_Remove_Patch),1) ;
    for ii = 1:length(P_Remove_Patch)
        D = C == P_Remove_Patch(ii) ;
        E = double(D).*AAA ;
        Patch_Mean_Biomass(ii) = mean(E,'all') ;
    end
    Poten_Loss_Biomass =  cumsum(Patch_Mean_Biomass) ;
    Difference = Poten_Loss_Biomass - Loss_Biomass ; 
    Check = find(Difference > 0) ;
    Remove_Patch = P_Remove_Patch(1:Check(1)) ;
    F = C ;
    F(ismember(C,Remove_Patch)) = 1000 ;
    G = F == 1000 ;
    H = 1- double(G) ;
    H(H == 0) = 1e-10 ;
    BBB = H.*AAA ;
end
