%% normalize for each metric
% for MMS: x = (x-min)/(max-min)
% input: MMS, is a matrix time × metric
MMS2=MMS
MMS2_min=zeros(1,7);
MMS2_max_min=zeros(1,7);
for i=1:7
    MMS2_min(1,i)=min(MMS2(:,i));
    MMS2_max_min(1,i)=max(MMS2(:,i))-min(MMS2(:,i));
    MMS2(:,i)=(MMS2(:,i)-MMS2_min(1,i))/MMS2_max_min(1,i);
end
%% normalize for each metric
% for SR: x = (x-min)/(max-min)
% input: SR, is a matrix time × metric
SR2=SR;
SR2_min=zeros(1,4);
SR2_max_min=zeros(1,4);
for i=1:4
    SR2_min(1,i)=min(SR2(:,i));
    SR2_max_min(1,i)=max(SR2(:,i))-min(SR2(:,i));
    SR2(:,i)=(SR2(:,i)-SR2_min(1,i))/SR2_max_min(1,i);
end
%% normalize for each metric
% for AQI: x = (x-min)/(max-min)
% input: AQI, is a matrix time × metric
AQI2=AQI;
AQI2_min=zeros(1,5);
AQI2_max_min=zeros(1,5);
for i=1:5
    AQI2_min(1,i)=min(AQI2(:,i));
    AQI2_max_min(1,i)=max(AQI2(:,i))-min(AQI2(:,i));
    AQI2(:,i)=(AQI2(:,i)-AQI2_min(1,i))/AQI2_max_min(1,i);
end

%% creat input tenor by a parameter T(cycle length)
% T: cycle length; W: window size
% performing only self-embedding transform for the first W slices
T=70;
index=1;
W=20;
index=1;
for i=1:(T*W-T+1)
%     M_mdt(:,index,:) = MMS2(i:i+T-1,:);
%     M_mdt(:,index,:) = SR2(i:i+T-1,:);
    M_mdt(:,index,:) = AQI2(i:i+T-1,:);
    index = index+1;
end
W_size=i
for i=T*W+1:T:5950
%     M_mdt(:,index,:)=MMS2(i:i+T-1,:);
%     M_mdt(:,index,:) = SR2(i:i+T-1,:);
    M_mdt(:,index,:) = AQI2(i:i+T-1,:);
    index=index+1;
end
