** using Stata to conduct the adjusted ROC analysis
** updated analysis after NG review 
** D Huo 7/1/2025

**************************************************************************************************************************************
 cd C:\Users\dhuo\Documents\tmp1
 global PRSvalid = "G:\BCAC/james.li/PRS_AABCG/output/ARISK_FAMH_SCORE_OUTPUT"
 global PRSvalidA = "G:\BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT"
 
** [1] Overall BC 
import delimited "$PRSvalidA/OVERALL/OVERALL_EUR_313_ALL_SNP.sscore", clear 
 ren score score_white
save OV_m0, replace 

import delimited "$PRSvalid/OVERALL_SINGLEANCESTRY_LDPRED2AUTO_hm3.sscore", clear 
 ren score score_single
save OV_m1, replace 

import delimited "$PRSvalid/OVERALL_XANCESTRY_PRScsx.sscore", clear 
 ren score score_x 
save OV_m2, replace 

import delimited "$PRSvalid/OVERALL_ENSEMBLE_GLMNET.sscore", clear 
 ren score score_e 
save OV_m3, replace 

import delimited "$PRSvalid/OVERALL_PROPTUNE_ALLSUBTYPE.sscore", clear 
 ren score score_p 
save OV_m4, replace 

import delimited "$PRSvalid/OVERALL_PROPTUNE_ALLSUBTYPE_WISDOM.sscore", clear 
 ren score score_p_wisdom 

merge 1:1 iid using OV_m0, gen(_m0) 
merge 1:1 iid using OV_m1, gen(_m1)
merge 1:1 iid using OV_m2, gen(_m2)
merge 1:1 iid using OV_m3, gen(_m3)
merge 1:1 iid using OV_m4, gen(_m4) 

pwcorr score_white score_single score_x score_e score_p score_p_wisdom
tabstat score_white score_single score_x score_e score_p score_p_wisdom, by(status) s(n mean sd)

 tab platform, gen(PLAT)

*Adjusted for covariate: single ancestry vs cross-ancestry models 
 rocreg status score_single score_x, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
   
*Adjusted for covariate: single ancestry vs Ensemble models 
 rocreg status score_single score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

*Adjusted for covariate: single ancestry vs best 3-subtype proportional-tuned models
 rocreg status score_single score_p, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

*Adjusted for covariate: best 3-subtype proportional-tuned vs. tuned Wisdom models
 rocreg status score_p_wisdom score_p, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
 
*Adjusted for covariate: cross-ancestry vs Ensemble models 
 rocreg status score_x score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
 
*Single ancestry model vs 313 model 
 rocreg status score_white score_single, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
* 313 model vs the best model 
 rocreg status score_white score_p, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
* 313 model vs the 2nd best model 
 rocreg status score_white score_p_wisdom, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
    
*Unadjusted
 roccomp status score_single score_x, summary 
 roccomp status score_single score_e, summary
 roccomp status score_x score_e, summary
 roccomp status score_white score_single score_x score_e score_p score_p_wisdom, summary 
 
 
** [2] ER+ BC 
import delimited "$PRSvalidA/ERPOS/ERPOS_EUR_313_ALL_SNP.sscore", clear 
 ren score score_white
save ERP_m0, replace 

import delimited "$PRSvalid/ERPOS_SINGLEANCESTRY_LDPRED2AUTO_hm3.sscore", clear 
 ren score score_single
save ERP_m1, replace 

import delimited "$PRSvalid/ERPOS_XANCESTRY_CTSLEB.sscore", clear 
 ren score score_x 
save ERP_m2, replace 

import delimited "$PRSvalid/ERPOS_ENSEMBLE_FSS.sscore", clear 
 ren score score_e 

merge 1:1 iid using ERP_m0, gen(_m0) 
merge 1:1 iid using ERP_m1, gen(_m1)
merge 1:1 iid using ERP_m2, gen(_m2)

pwcorr score_white score_single score_x score_e  
tabstat score_white score_single score_x score_e, by(status) s(n mean sd)

quietly: tab platform, gen(PLAT)

*Adjusted for covariate: single ancestry vs cross-ancestry models 
 rocreg status score_single score_x, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
   
*Adjusted for covariate: single ancestry vs Ensemble models 
 rocreg status score_single score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

*Adjusted for covariate: cross-ancestry vs Ensemble models 
 rocreg status score_x score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

* 313-SNP model vs single ancestry
 rocreg status score_white score_single, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
   
*313-SNP model vs Ensemble models 
 rocreg status score_white score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

*Unadjusted
 roccomp status score_single score_x, summary 
 roccomp status score_single score_e, summary
 roccomp status score_x score_e, summary
 roccomp status score_white score_single score_x score_e, summary 
 
  
** [3] ER- BC 
import delimited "$PRSvalidA/ERNEG/ERNEG_EUR_313_ALL_SNP.sscore", clear 
 ren score score_white
save ERN_m0, replace 
 
import delimited "$PRSvalid/ERNEG_SINGLEANCESTRY_LDPRED2AUTO_hm3.sscore", clear 
 ren score score_single
save ERN_m1, replace 

import delimited "$PRSvalid/ERNEG_XANCESTRY_PLINKCT_0.1.sscore", clear 
 ren score score_x 
save ERN_m2, replace 

import delimited "$PRSvalid/ERNEG_ENSEMBLE_GLMNET.sscore", clear 
 ren score score_e 
merge 1:1 iid using ERN_m0, gen(_m0)
merge 1:1 iid using ERN_m1, gen(_m1)
merge 1:1 iid using ERN_m2, gen(_m2)

pwcorr score_white score_single score_x score_e  
tabstat score_white score_single score_x score_e, by(status) s(n mean sd)

quietly: tab platform, gen(PLAT)

*Adjusted for covariate: single ancestry vs cross-ancestry models 
 rocreg status score_single score_x, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
   
*Adjusted for covariate: single ancestry vs Ensemble models 
 rocreg status score_single score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

*Adjusted for covariate: cross-ancestry vs Ensemble models 
 rocreg status score_x score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

* 313-SNP model vs single ancestry
 rocreg status score_white score_single, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
   
*313-SNP model vs Ensemble models 
 rocreg status score_white score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

 
  
** [4] TNBC 
import delimited "$PRSvalidA/TNBC/TNBC_EUR_330_ALL_SNP.sscore", clear 
 ren score score_white
save TN_m0, replace 

import delimited "$PRSvalid/TNBC_SINGLEANCESTRY_LASSOSUM2_hm3.sscore", clear 
 ren score score_single
save TN_m1, replace 

** cannot find the Plink C+T all variant model, so choose another model with similar performanc (AUC=0.639)
import delimited "$PRSvalid/TNBC_XANCESTRY_PLINKCT_0.2.sscore", clear   
 ren score score_x 
save TN_m2, replace 

** cannot find the Plink C+T all variant model, assume it is the PRSICE2 model (AUC=0.639)
import delimited "$PRSvalid/TNBC_XANCESTRY_PRSICE2.sscore", clear   
 ren score score_x2 
save TN_m2A, replace 

import delimited "$PRSvalid/TNBC_ENSEMBLE_FSS.sscore", clear 
 ren score score_e_FSS  
save TN_m3, replace 

import delimited "$PRSvalid/TNBC_ENSEMBLE_GLMNET.sscore", clear 
 ren score score_e_GLM 

merge 1:1 iid using TN_m0, gen(_m0) 
merge 1:1 iid using TN_m1, gen(_m1)
merge 1:1 iid using TN_m2, gen(_m2)
merge 1:1 iid using TN_m2A, gen(_m2A)
merge 1:1 iid using TN_m3, gen(_m3)

quietly: tab platform, gen(PLAT)

pwcorr score*
drop score_e_GLM score_x 
tabstat score_*, by(status) s(n mean sd)

 *Adjusted for covariate: single ancestry vs cross-ancestry Price-2 model
 rocreg status score_single score_x2, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
 
*Adjusted for covariate: single ancestry vs Ensemble models 
 rocreg status score_single score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31433)	
 * rocreg status score_single score_e, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) nodots breps(1000) bseed(31433)	

*Adjusted for covariate: cross-ancestry vs Ensemble models 
 rocreg status score_x score_e_, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
 
* 313-SNP model vs single ancestry
 rocreg status score_white score_single, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	
   
* 313-SNP model vs cross-ancestry Price-2 model 
 rocreg status score_white score_x2, ctrlcov(PLAT1-PLAT6 age pc1-pc10) ctrlmodel(linear) pvc(normal) nodots breps(1000) bseed(31415923)	

 
*Unadjusted
 roccomp status score_single score_x, summary 
 roccomp status score_single score_e, summary
 roccomp status score_x score_e, summary
 roccomp status score_*, summary