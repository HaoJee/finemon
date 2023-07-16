function [W,U_W_t,U_W_index,ischange_U_W,rank_W]=paper_subroutine_sw_removeFromW(W,U_W_t,U_W_index,count_map,ischange_U_W,rank_W,yita,i)
    oldest_slice_index=W(1);
    W=W(1,2:length(W));
    if ismember(oldest_slice_index,U_W_index)
        if isKey(count_map,oldest_slice_index)
            count=count_map(oldest_slice_index);
            U_W_remove=[];
            if size(U_W_t,2)>count
                U_W_remove=U_W_t(:,1:count,:);
                U_W_t=U_W_t(:,count+1:size(U_W_t,2),:);
                
                pesudoinverse_U_W=tpinv(U_W_t);
                Pu_slideW2=tprod(U_W_t,pesudoinverse_U_W);
                for j=1:count
                    proRes=tprod(Pu_slideW2,U_W_remove(:,j,:));
                    estimator=(norm(tensor(U_W_remove(:,j,:)-proRes))^2)/(norm(tensor(U_W_remove(:,j,:)))^2);
                    if estimator>yita
                        ischange_U_W=1;
                        rank_W=rank_W-1;
                    end
                end
                U_W_index=setdiff(U_W_index,oldest_slice_index);
            end
        else
            U_W_remove=[];
            if size(U_W_t,2)>1
                U_W_remove(:,1,:)=U_W_t(:,1,:);
                U_W_t=U_W_t(:,2:size(U_W_t,2),:);
%                 ischange_U_W=1;
                proRes=tprod(tprod(U_W_t,tpinv(U_W_t)),U_W_remove);
                estimator=(norm(tensor(U_W_remove-proRes))^2)/(norm(tensor(U_W_remove))^2);
                if estimator>yita
                    ischange_U_W=1;
                    rank_W=rank_W-1;
                end
                U_W_index=setdiff(U_W_index,oldest_slice_index);
            end
        end
    end

%     窗口移动，更新pareto参数
%     for pd_i = 1:K
%         M_pd_i(:,:) = M(:, W, pd_i);
%         pareto_delta(pd_i) = mean(M_pd_i(:));
%         pareto_x_min(pd_i) = min(M_pd_i(:));
%     end
    rank_W=rank_W+1;
    W=[W i];
end