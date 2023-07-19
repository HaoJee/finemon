# FineMon
### Introduction：

This is a description of our monitoring scheme-FineMon. The code of FineMon is the MATLAB version. 

The coding process is shown below figure:

![finemon_process](finemon_process.png)

We introduce each code files as follows:

data_proprocess.m:  Data preprocessing, where the original data is processed into an input tensor according to the flow in the Figure.

FineMon.m:  Main file, this file can perform FineMon for each dataset (input data) and calculate the evaluation criteria (NMAE, Cos, Sample Ratio).

getPerformanceNC_orign.m: Calculate all the data's NMAE and Cos.

subroutine_residual.m: Calculate the residual(.) mensioned in equation (3).

subroutine_sw_removeFromW.m: Update the active window by removing the oldest slice, calculate the common window's rank, and estimate the active window's rank by the upper bound in the Theorem 4.1.

subroutine_sw_STupdateU.m: Perform self-embedding transform for full sampling slices, this is proposed in Section 4.3.1

subroutine_delayRecovery.m: The delay-recovery strategy mensioned in Section 4.3.

subroutine_delayRecovery_all.m: Perform the delay-recovery strategy for all previous incompleted slices.





library function:

teye.m:  an identity tensor

tinv.m: the inverse of tensor

tpinv.m: the pseudo inverse of tensor

tprod.m: the tensor product "*"

tran.m: the transpose of tensor

tsvd.m: the t-SVD decomposition





Dataset:

MMS: A server performance dataset, the original link is [QAZASDEDC/TopoMAD: A Spatio-Temporal Deep Learning Approach for Unsupervised Anomaly Detection in Cloud Systems (TNNLS) (github.com)](https://github.com/QAZASDEDC/TopoMAD).   The original data is OriginalData_MMS_mat. The input tensor processed by data_preprocess.m is InputTensor_MMS_T=40_W=8.mat.

AQI:  the original link is https://www.microsoft.com/en-us/research/uploads/prod/2016/02/Data-1.zip  The original data is OriginalData_AQI(precompletion_by_LRTC)_mat. The input tensor processed by data_preprocess.m is InputTensor_AQI_T=50_W=15.mat.

SR: Our traffic statistics telemetry data (company name withheld to comply with anonymity rules). The original data is OriginalData_SR.mat. The input tensor processed by data_preprocess.m is InputTensor_SR_T=50_W=15.mat.

