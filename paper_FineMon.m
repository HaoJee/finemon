%ʹ�ò����U
%ʹ�� U union
% W_size=991;%MMS-23-45
W_size=281;

% 15  0.03 1e-21
repeatTimes=1;
theta=2.5e-3;%MMS
yita=4e-29;
% theta=2.4e-4;%SR
% yita=6e-28;
% theta=0.5e-3;%AQI
% yita=5e-29;
% theta=9e-4;%AQI5
% yita=4e-29;  
alpha=2;

I=size(M,1);
J=size(M,2);
K=size(M,3);

            [R,omega,h_incoms,estimators,rs,ms,com,U_W_t,U_union_t,U_W_index,times,times_update_rs,times_estimatorFull,times_completion,times_compute_m,times_oumega,times_fullcompletion,times_fullfor,OptimizationIDX, pareto_oemga]=tensorCompletion_sw(M,W_size,I,J,K,theta,yita,alpha);
            
            %X_res=X_res+X;
            %incoms_repeatRes=[incoms_repeatRes incoms];
            R(find(R<0))=0;
            Rm=[];
            for j=W_size+1:J
                Rm(size(Rm,1)+1:size(Rm,1)+I,:)=R(:,j,:);
            end
        
        %  %%*******case2-1 R��M ���ǰ�ָ���źŹ�һ�������ݣ�������һ��
%             %     MMS2_max_min  MMS2_min
%             %�ȷ���һ��
%             %orig_data_M=MMS(751:3600,:);%MMS
%             orig_data_M=[];
%             for j=W_size+1:J
%                 orig_data_M(size(orig_data_M,1)+1:size(orig_data_M,1)+I,:)=M(:,j,:);
%             end
%             p=0.2;
%             [p_NMAEs,p_COSes,p_peakCoverages]=getPerformanceNCP_orign(orig_data_M, Rm, p);
%             
%             %�ټ���ָ��
%             %ʹ��δ����λ��
%             unOmega=ones(I,J,K)-omega;
%             p_unOmega_NMAEs=zeros(1,K);
%             for k=1:K
%                 Mk = orig_data_M(:,k);
%                 RMk = Rm(:,k);
%                 unOmegak(:,:) = unOmega(:,W_size+1:J,k);
%                 fenmuk = sum(abs(Mk.*unOmegak(:)))
%                 if fenmuk>0
%                     p_unOmega_NMAEs(1,k) = sum(abs(RMk.*unOmegak(:)-Mk(:).*unOmegak(:)))/fenmuk;
%                 else
%                     p_unOmega_NMAEs(1,k) = 0;
%                 end
%             end

             %  %%*******case2 R��M ���ǰ�ָ���źŹ�һ�������ݣ���Ҫ����һ����ȥ
            %     MMS2_max_min  MMS2_min
            %�ȷ���һ��
            N_max_min=MMS2_max_min;
            N_min=MMS2_min;
%             N_max_min=KPIs2_max_min;
%             N_min=KPIs2_min;
%             N_max_min=Mone2_max_min;
%             N_min=Mone2_min;
            RM_orign=[];
            for k=1:K
                RM_orign(:,k) = Rm(:,k)*N_max_min(k)+N_min(k);
            end
            orig_data_M=MMS(8*40+1:(J-W_size+8)*40,:);%MMS
%             orig_data_M=KPIs(15*50+1:(J-W_size+15)*50,:);%SR
%             orig_data_M=Mone(15*50+1:(J-W_size+15)*50,:);%AQI
            p=0.1;
            [p_NMAEs,p_COSes,p_peakCoverages]=getPerformanceNCP_orign(orig_data_M, RM_orign, p);
            
            %�ټ���ָ��
            %ʹ��δ����λ�õĹ�һ������
            unOmega=ones(I,J,K)-omega;
            p_unOmega_NMAEs=zeros(1,K);
            for k=1:K
                Mk = orig_data_M(:,k);
                RMk = RM_orign(:,k);
                unOmegak(:,:) = unOmega(:,W_size+1:J,k);
                fenmuk = sum(abs(Mk.*unOmegak(:)));
                if fenmuk>0
                    p_unOmega_NMAEs(1,k) = sum(abs(RMk.*unOmegak(:)-Mk(:).*unOmegak(:)))/fenmuk;
                else
                    p_unOmega_NMAEs(1,k) = 0;
                end
            end
        
            %�ܵ�NMAE
            NMAE=sum(sum(abs(orig_data_M-Rm)))/sum(sum(abs(orig_data_M)));
            unIDX=[];
            for j=W_size+1:J
                unIDX(size(unIDX,1)+1:size(unIDX,1)+I,:)=unOmega(:,j,:)
            end
            NMAE_unomega=sum(sum(abs((orig_data_M-Rm).*unIDX)))/sum(sum(abs(orig_data_M.*unIDX)));
            
            
            %�ܲ�����
            samples_total=sum(sum(sum(omega)));
            sampleRatio_total=samples_total/(I*J*K);
            
            %�������ڽ׶β�����
            samples_sw=samples_total-W_size*I*K;
            %     samples_sw=samples_total-W_size*I*J;
            %     sampleRatio_sw=sum(sum(sum(omega(:,:,W_size+1:K))))/(I*J*(K-W_size));
            sampleRatio_sw=samples_sw/(I*(J-W_size)*K);
            
            %ÿһ��Ĳ����������ȼ�����������Ա�ʾƵ��
            perSliceSamples=zeros(1,K);
            perSliceSamples_tube=sum(sum(omega));
            for j=1:K
                perSliceSamples(1,j)=perSliceSamples_tube(:,:,j);
            end

function [R,omega,h_incoms,estimators,rs,ms,com,U_W_t,U_union_t,U_W_index,times,times_update_rs,times_estimatorFull,times_completion,times_compute_m,times_oumega,times_fullcompletion,times_fullfor,OptimizationIDX, pareto_oemga]=tensorCompletion_sw(M,W_size,I,J,K,theta,yita,alpha)
    %��ʼ��
    W=[];
    U_W_t=[];
    U_W_t_indep=[];%���Զ���Uw�����ټ��㿪��
    U_union_t=[];
    U_W_index=[];
    U_W_t_indep_index=[];
    rank_W=0;
    u1=1;

    com=zeros(1,J);
    count_map=containers.Map;
    estimators=zeros(1,J);

    incomplete_spaceSlice=0;
    h_incoms=[];

    R=[];
    omega=[];
    
    Pu=[];

    pesudoinverse_U_W_omega_t=[];
    Pu_omega=[];

    rs=zeros(1,J);
    ms=zeros(1,J);
    times=zeros(1,J);%��ʱ�� (ÿ�����ݴӲ���������ʱ�䣬���ܶ�����)
    times_update_rs=zeros(1,J);%�������ڣ�����r��ʱ��
    times_estimatorFull=zeros(1,J);%����ȫ������ʱ��
    times_completion=zeros(1,J);%���ֲ���������ʱ��
    times_compute_m=zeros(1,J);
    times_oumega=zeros(1,J);

    times_fullcompletion=zeros(1,J);%ȫ�����׶Σ����һ���ʱ��
    times_fullfor=zeros(1,J);%ȫ�����׶Σ�ִ��forѭ����ʱ��

    times_train2=zeros(1,J);

    incomplete_spaceSlices=zeros(1,J);

    alpha_count=0;
    Mx=[];

    ischange_U_W=1;%�ж�U�Ƿ�仯�����û�䣨��û�䣩�����������ʱ�Ϳ���ʡȥ���㿪��
    isSame_UWTindep=0;
    
    OptimizationIDX=zeros(1,J);%ͳ��ʱ���Ż����ĸ�slice������������㡣

%     TVS������ ��ֵ����\delta, ��Сֵ����x_m, ����alpha_m, L
    pareto_delta = zeros(K,1);
    pareto_x_min = zeros(K,1);
    pareto_alpha = zeros(K,1);
    pareto_L = 0;
    pareto_oemga=zeros(I,J,K);

    for i=1:J
        if length(W)<W_size
            %ѵ���׶Σ�ȫ����
            st_train=tic;
            W=[W i];
            Mi=M(:,i,:);
            R(:,i,:)=Mi;
            com(1,i)=-1;
            omega(:,i,:)=ones(I,K);
            ms(1,i)=I*K;

            if size(U_W_t_indep,2)==0
                estimator=1;
            else
                %Ϊ�˼���U_W_t�Ĺ�ģ��ʹ�����Զ�����U_W_t_indep
                st_train2=tic;
                if ~isSame_UWTindep
                    Pu_train2=tprod(U_W_t_indep,tpinv(U_W_t_indep))
                end
                estimator=(norm(tensor(Mi-tprod(Pu_train2,Mi)))^2)/(norm(tensor(Mi))^2);
                times_train2(1,i)=toc(st_train2);
            end
            normalize=norm(Mi(:));
            if normalize==0
                tempU=Mi;
            else
                tempU=Mi/normalize;
            end
            if size(U_W_t,2)==0
                n=1;
            else
                n=size(U_W_t,2)+1;
            end
            U_W_t(:,n,:)=tempU;
            U_W_index=[U_W_index i];
            estimators(1,i)=estimator;
            if estimator>yita
                U_W_t_indep(:,size(U_W_t_indep,2)+1,:)=tempU;%���Զ�����Uw
                U_W_t_indep_index=[U_W_t_indep_index i];
                isSame_UWTindep=0;
                rank_W=rank_W+1;
            else
                isSame_UWTindep=1;%���Ա����ظ�����Pu_train2
            end
            time=toc(st_train);
            times(1,i)=time;
            rs(1,i)=rank_W;
        else
            %�������ڽ׶�
            %ʱ�䣺times_r��times_UWtr+����ͶӰ   vs    times_r������ͶӰ
            st_slideW=tic;
            %subroutine_sw_removeFromW���������ڣ���U���Ƴ����ϵ��У�����r������r'
            %------�������������W��U_W,U_W_t,U_W_index,count_map,ischange_U_W,rank_W
            %------������º�Ĳ�����W,U_W,U_W_t,U_W_index,ischange_U_W,rank_W
%             [W,U_W_t,U_W_index,ischange_U_W,rank_W]=subroutine_sw_removeFromW(W,U_W_t,U_W_index,count_map,ischange_U_W,rank_W,yita);
            [W,U_W_t,U_W_index,ischange_U_W,rank_W]=paper_subroutine_sw_removeFromW(W,U_W_t,U_W_index,count_map,ischange_U_W,rank_W,yita,i);
            times_r=toc(st_slideW);
            times_update_rs(1,i)=times_r;

            if incomplete_spaceSlice>0
                st_full1=tic;
                Mi=M(:,i,:);
                R(:,i,:)=Mi;
                omega(:,i,:)=ones(I,K);
                ms(1,i)=I*K;
                com(:,i)=-1;
                alpha_count=alpha_count+1;
                Mx(size(Mx,1)+1:size(Mx,1)+I,:,:)=Mi;
                if alpha_count==alpha
                    alpha_count=0;
                    %subroutine_sw_MDTupdateU��MDT����ȫ�������ݣ�����Uw��Uall
                    %-----���������Mx(MDT����������),U_W,U_W_t,U_all,U_all_t,ischange_U_W,estimators
                    %-----������º��U������������U_W,U_W_t,U_all,U_all_t,ischange_U_W,estimators,count
%                     [U_W_t,U_all_t,ischange_U_W,estimators,count]=subroutine_sw_MDTupdateU(Mx,U_W_t,U_all_t,i,I,ischange_U_W,estimators,yita);
                    [U_W_t,U_union_t,ischange_U_W,estimators,count]=paper_subroutine_sw_MDTupdateU(Mx,U_W_t,U_union_t,i,I,ischange_U_W,estimators,yita);
                    
                    time=toc(st_full1);
                    times(1,i)=time;
                    if count>0
                        rank_W=rank_W+count;
                        U_W_index=[U_W_index i];
                        keys=count_map.keys;
                        keys=[keys i];
                        values=count_map.values;
                        values=[values count];
                        count_map=containers.Map(keys,values);
                    end
                    Mx=[];

                    %subroutine_tensorcompletion_full������incomplete_spaceSlice��
                    %------���������M,R,U_all,U_all_t,omega,h_incoms,com,incomplete_spaceSlice,times,times_fullcompletion
                    %------�������������Լ����º�ı�����R,com,h_incoms,incomplete_spaceSlice,times,times_fullcompletion
%                     [R,com,h_incoms,incomplete_spaceSlice,times,times_fullcompletion]=subroutine_tensorcompletion_full(M,R,U_all_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullcompletion,yita);
                    [R,com,h_incoms,incomplete_spaceSlice,times,times_fullcompletion]=paper_subroutine_tensorcompletion_full(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullcompletion,yita);
           
                    %subroutine_tensorcompletion_fullfor������h_incoms���������
                    %------���������M,R,U_all,U_all_t,omega,h_incoms,com,incomplete_spaceSlice,times,times_fullfor
                    %------�������������Լ����±�����R,com,h_incoms,times,times_fullfor
%                     [R,com,h_incoms,times,times_fullfor]=subroutine_tensorcompletion_fullfor(M,R,U_all_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullfor,yita);
                    [R,com,h_incoms,times,times_fullfor]=paper_subroutine_tensorcompletion_fullfor(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullfor,yita);
                    incomplete_spaceSlice=0;

                    if length(h_incoms)==0
                        U_union_t=[];
                    end
                else
                    time=toc(st_full1);
                    times(1,i)=time;
                end

            else%if incomplete_spaceSlice=0
                st_sample=tic;
                %���в��ֲ���
                if rank_W==1
                    m=min(ceil(theta*rank_W*rank_W),I*K);
                else
                    m=min(ceil(theta*rank_W*rank_W*log(rank_W)),I*K);
                end
                times_compute_m(1,i)=toc(st_sample);
                m=floor(m/K);%ƽ����ÿһ��
                ms(1,i)=m*K;
                %������Ĳ�������
                omega_column=Get_Array_equalInterval(I,m);
                %ÿһ�еĲ���λ������
                omega_slice_index=find(omega_column==1);
                %����ֵ
                if m~=I
                    for j=1:length(omega_slice_index)
                        if omega_slice_index(j)<I
                            omega(omega_slice_index(j),i,:)=1;
                            for k = 1:K
                                MWk = M(:,W,k);
                                omegaWk = omega(:,W,k);
                                mWk = MWk(:);
                                oWk = omegaWk(:);
                                oWk_IDX = find(oWk == 1);
                                pareto_delta(k)=3*mean(mWk(oWk_IDX));
                                mWk_oWk = mWk(oWk_IDX);
                                pareto_x_min(k)=min(mWk_oWk(mWk_oWk>0));
                                pareto_alpha(k) = length(oWk_IDX)/(sum( log( mWk(oWk_IDX)/pareto_x_min(k) ) ));
                            end
                            m_ji = M(omega_slice_index(j),i,:);
                            is_large = m_ji > pareto_delta;
                            if max(is_large)==1
                                % ��⵽��ֵ������L
                                min_pareto_delta = min(pareto_delta);
                                pareto_x_min_delta = pareto_x_min/min_pareto_delta;
                                pareto_x_min_delta_alpha = power(pareto_x_min_delta,pareto_alpha);
                                L =  floor( min(I-omega_slice_index(j), (ceil(I/m)-1))*max(pareto_x_min_delta_alpha) );
    
                                %��(omega_slice_index(j), omega_slice_index(j+1) or I) �ɼ� L������
                                Ll = min(I-omega_slice_index(j), (ceil(I/m)-1));
                                if Ll<L
                                    L=Ll;
                                end
                                omega_column_L = Get_Array_equalInterval(Ll,L);
                                omega_column( omega_slice_index(j)+1 : omega_slice_index(j)+Ll ) = omega_column_L;
                                omega( omega_slice_index(j) + find(omega_column_L==1), i, :) = 1;
                                pareto_oemga(omega_slice_index(j)+find(omega_column_L==1),i,:) = 1;
%                                 if length(omega_slice_index)>j
%                                     omega_column(omega_slice_index(j)+1: omega_slice_index(j+1)-1) = 1;
%                                     omega(omega_slice_index(j)+1: omega_slice_index(j+1)-1 ,i,:)=1;
%                                 else
%                                     omega_column(omega_slice_index(j)+1: I) = 1;
%                                     omega(omega_slice_index(j)+1: I ,i,:)=1;
%                                 end
                            end
                        end
                    end
                end
                for j=1:K
                    omega_slice(:,j)=omega_column;
                end
                omega(:,i,:)=omega_slice;
                 %ÿһ�еĲ���λ������
                omega_slice_index=find(omega_column==1);
                Mi_omega_1(:,:)=M(:,i,:);
                Mi_omega_1=Mi_omega_1.*omega_slice;
                Mi_omega=[];
                Mi_omega(:,1,:)=Mi_omega_1;
                times_oumega(1,i)=toc(st_sample);
                
%                 Mi_omega_t=[];
%                 Mi_omega_t(:,1,:)=Mi_omega;
                R(:,i,:)=Mi_omega;
                st_estimatorFull=tic;
                %����estimator
                %������λ����ͬ��Uͬ�ϴ�û�б仯������Ҫ�ٴμ���
                isSame_sample=0;
                %subroutine_systemDecision_isCompletable����ⲿ�ֲ����Ƿ�����
                %-----���������Mi_omega(���ֲ�������),Mi_omega_t,omega,omega_slice,i,U_W(�ӿռ�),pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W
                %-----���estimator�Լ����º�ı�����estimator,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample
%                 [estimator,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample]=subroutine_systemDecision_isCompletable(Mi_omega_t,omega,omega_slice,i,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W,W_size);
                [estimator,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,OptimizationIDX]=paper_subroutine_systemDecision_isCompletable(Mi_omega,omega,omega_slice,omega_slice_index,i,U_W_t,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W,W_size,OptimizationIDX);

                times_estimatorFull(1,i)=toc(st_estimatorFull);
                estimators(1,i)=estimator;
                %���
                if estimator>yita
                    st_estimatorFull2=tic;
                    if m~=I
                        incomplete_spaceSlice=i;
                        h_incoms=[h_incoms i];
                        incomplete_omega_column=omega_column;
                        alpha_count=0;
    
                        if size(U_union_t,2)==0
                            U_union_t = U_W_t;
                        end
                    else
                        %����U��U_union:  ����ֵ==��ʵ��
                        normalize=norm(Mi_omega_1(:));
                        if normalize==0
                            tempU=Mi_omega_1;
                        else
                            tempU=Mi_omega_1./normalize;
                        end
                        U_W_t(:,size(U_W_t,2)+1,:)=tempU;
                        U_W_index=[U_W_index i];
                        if size(U_union_t,2)>0
                            U_union_t(:,size(U_union_t,2)+1,:)=tempU;
                        end
                    end

                    times_estimatorFull(1,i)=times_estimatorFull(1,i)+toc(st_estimatorFull2);
                else
                    st_completion=tic;

                    if i==W_size+1 || ~isSame_sample
                        %ridge regression
                        U_W_omega=U_W_t(omega_slice_index,:,:);
                        ridge = tprod(tran(U_W_omega),U_W_omega);
                        ridge_eye = 0.01*teye(size(ridge,1), size(ridge,3));
                        pesudoinverse_U_W_omega_t=tprod(tinv(ridge-ridge_eye),tran(U_W_omega)); 

                        Pu=tprod(U_W_t,pesudoinverse_U_W_omega_t);
%                         ischange_U_W=0;%ÿ�μ���Uû�䣬�����ˣ���ֵ����Ϊ1
                        ischange_U_W=1;%ÿ��U��
                    else
                        OptimizationIDX(1,i)=OptimizationIDX(1,i)+1;
                    end
                    %subroutine_tensorcompletion��ͨ��Uw��䲿�ֲ�������
                    %------���������Pu(Uw),Mi_omega(���ֲ�������),i,omega_slice_index
                    %------������������R�͸��±���com
%                     [R,com]=subroutine_tensorcompletion(R,Pu,Mi_omega_t,Mi_omega,i,omega_slice_index);
                    [R,com]=paper_subroutine_tensorcompletion(R,Pu,Mi_omega,i,omega_slice_index,com);
                    times_completion(1,i)=toc(st_completion);
                    rank_W=rank_W-1;
                end

                time=toc(st_sample);
                times(1,i)=time;
            end %if incomplete_spaceSlice>0
            rs(1,i)=rank_W;
        end
        incomplete_spaceSlices(1,i)=incomplete_spaceSlice;

    end%end for
end