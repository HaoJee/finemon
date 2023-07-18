# finemon
This is a descripution about our monitoring scheme-FineMon. The code of FineMon is MATLAB version. 

The code process is shown as below figure:

![finemon_process](E:\workplace\matlab\multi-metric\current\FineMon\finemon_process.png)

We introduce each code files as follows:

data_proprocess.m:  Data preprocessing, where the original data is processed into an input tensor according to the flow in the Figure.

FineMon.m:  Main file, this file can perform FineMon for each dataset (input data) and calculate the evalution criteria (NMAE, Cos, Sample Ratio).

getPerformanceNC_orign.m: Calculate the all data's NMAE and Cos.

subroutine_residual.m: Calculate the residual(.) mensioned in equation (3).

subroutine_sw_removeFromW.m: Update the active window by removing the oldest slice, calculate the common window's rank, and estimate the active window's rank by the upper bound in the Theorem 4.1.

subroutine_sw_STupdateU.m: Perform self-embedding transform for full sampling slices, this is proposed in Section 4.3.1

subroutine_delayRecovery.m: The delay-recovery strategy mensioned  in Section 4.3.

subroutine_delayRecovery_all.m: Perform the delay-recovery strategy for all previous incompleted slices.





library function:

teye.m:  an identity tensor

tinv.m: the inverse of tensor

tpinv.m: the pesudo inverse of tensor

tprod.m: the tensor product "*"

tran.m: the transpose of tensor

tsvd.m: the t-SVD decomposition





