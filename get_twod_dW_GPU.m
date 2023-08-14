function [dW1,dW2] = get_twod_dW_GPU(bj,kappa)
    J = size(bj) ;
    if kappa == 1 % get xi_j
        nnr = randn(J(1),J(2),'gpuArray') ; 
        nnc = randn(J(1),J(2),'gpuArray') ;
    else % sum over kappa steps
        nnr = squeeze(sum(randn(J(1),J(2),kappa,'gpuArray'),4)) ;
        nnc = squeeze(sum(randn(J(1),J(2),kappa,'gpuArray'),4)) ;
    end
    nn2 = nnr + 1i*nnc ; tmphat = bj.*nn2 ;
    tmp = ifft2(tmphat) ; dW1 = real(tmp) ; dW2 = imag(tmp) ;
end