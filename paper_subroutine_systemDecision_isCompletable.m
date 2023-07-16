function [estimator,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,OptimizationIDX]=paper_subroutine_systemDecision_isCompletable(Mi_omega,omega,omega_slice,omega_slice_index,i,U_W,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W,W_size,OptimizationIDX)
    if size(U_W,2)==0
        estimator=1;
    else
        %若采样位置相同且U同上次没有变化，则不需要计算
        isSameOmega=0;
        for j=1:size(omega_slice,1)
            if omega_slice(j,1)==omega(j,i-1,1)
                isSameOmega=1;
            else
                isSameOmega=0;
                break
            end
        end
        isSame_sample=(~ischange_U_W) && isSameOmega;
        isSame_sample=0;
        if i==W_size+1 || ~isSame_sample
            if length(omega_slice_index)>0
                U_W_omega=U_W(omega_slice_index,:,:);
                pesudoinverse_U_W_omega_t=tpinv(U_W_omega);
%                 pesudoinverse_U_W_omega_t=tprod(tinv(tprod(tran(U_W_omega),U_W_omega)),tran(U_W_omega)); 
%                 ridge = tprod(tran(U_W_omega),U_W_omega);
%                 ridge_eye = teye(size(ridge,1), size(ridge,3));
%                 pesudoinverse_U_W_omega_t=tprod(tinv(ridge-ridge_eye),tran(U_W_omega)); 
                Pu_omega=tprod(U_W_omega,pesudoinverse_U_W_omega_t);
            end
        else
            OptimizationIDX(1,i)=OptimizationIDX(1,i)+1;
        end
        proRes=[];
        if length(omega_slice_index)>0
            Mi_omega_1=Mi_omega(omega_slice_index,:,:);
            proRes(:,:)=tprod(Pu_omega,Mi_omega_1);
            estimator=(norm(Mi_omega_1(:,:)-proRes,'fro')^2)/(norm(Mi_omega_1(:,:),'fro')^2);
        else
            estimator=0;
        end
    end
end