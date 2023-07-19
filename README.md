# finemon
<<<<<<< HEAD
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





=======
code describe：
file:
   1. paper_FineMon.m  is the main code for our paper.
   2. paper_subrotine_sw_MDTupdateU.m:   a subfunction in FineMon,  update the subspace in the "delay-recovery strategy" in our paper
   3. paper_subrotine_sw_removeFromW.m:   a subfunction in FineMon,  update the subspace by removing the first slice
   4. paper_subrotine_systemDecision_isCompletable.m:   a subfunction in FineMon,  complete the residual(.) to determine if the sampled slice can be recovery
   5. paper_subrotine_tensorcompletion.m:   a subfunction in FineMon,  recover the un-sampled data from sampled slice by  ESTC
   6. teye.m : base function for tensor,   create an identity tensor
   7. tinv.m: base function for tensor,  compute tensor inverse
   8. tpinv.m: base function for tensor, compute the pesudo-inverse of tensor
   9. tprod.m: base function for tensor, "*" in paper, compute the tensor-tensor product,
   10. tran.m: base function for tensor, compute the translation of tensor
   11. tsvd.m base function for tensor, compute the t-SVD decomposition of a tensor
>>>>>>> 72513ef30d1442ae3b8547de414fc96d6c6150d3
