# finemon
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
