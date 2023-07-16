%Mi_omega：表示包含0的采样数据
function [R,com]=paper_subroutine_tensorcompletion(R,Pu,Mi_omega,i,omega_slice_index,com)
    Mi_omega_1=Mi_omega(omega_slice_index,1,:);
    R(:,i,:)=tprod(Pu,Mi_omega_1);
    R(omega_slice_index,i,:)=Mi_omega(omega_slice_index,1,:);
    com(:,i)=1;
end