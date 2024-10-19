Supplementary information / reproducible research files for the manuscript

1. Introduction
----------------------------------------------------------------------------------------------------- 
Title: "A Two-stage Screening and Selection Approach for Ultra-High Dimensional Competing Risks Data"

Authors: Li,Erqian, Xiong,Wei, Pan,Han and Tian,Maozai
Authors of the code: Li,Erqian, Xiong,Wei(xiongwei@uibe.edu.cn), Pan,Han

2. Configurations
-----------------------------------------------------------------------------------------------------
The code uses nine additional R packages: survival, mvtnorm, lars, crrp, MFSIS, cmprsk, riskRegression, pec, evd.
These packages can be installed using the devtools package and the following commands:

install.packages("survival")
install.packages("mvtnorm")
install.packages("lars")
install.packages("crrp")
install.packages("MFSIS")
install.packages("cmprsk")
install.packages("riskRegression")
install.packages("pec")
install.packages("evd")

The code was written/evaluated in RStudio with the following software versions:

R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pec_2023.04.12            prodlim_2023.08.28        riskRegression_2023.12.21
 [4] MFSIS_0.2.0               crrp_1.0                  cmprsk_2.2-11            
 [7] Matrix_1.6-1.1            MSBVAR_0.9-3              lars_1.3                 
[10] mvtnorm_1.2-2             survival_3.5-5           

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0    dplyr_1.1.2         fastmap_1.1.1      
 [4] TH.data_1.1-2       digest_0.6.33       rpart_4.1.19       
 [7] lifecycle_1.0.4     cluster_2.1.4       magrittr_2.0.3     
[10] compiler_4.3.1      rlang_1.1.1         Hmisc_5.1-0        
[13] tools_4.3.1         utf8_1.2.3          data.table_1.14.8  
[16] knitr_1.43          timereg_2.0.5       htmlwidgets_1.6.2  
[19] bit_4.0.5           reticulate_1.30     multcomp_1.4-25    
[22] KernSmooth_2.23-21  polspline_1.1.24    foreign_0.8-84     
[25] numDeriv_2016.8-1.1 Ball_1.3.13         nnet_7.3-19        
[28] grid_4.3.1          fansi_1.0.4         mets_1.3.3         
[31] xtable_1.8-4        colorspace_2.1-0    future_1.33.0      
[34] ggplot2_3.4.3       globals_0.16.2      scales_1.2.1       
[37] iterators_1.0.14    MASS_7.3-60         cli_3.6.1          
[40] rmarkdown_2.23      crayon_1.5.2        rms_6.7-1          
[43] generics_0.1.3      rstudioapi_0.15.0   future.apply_1.11.0
[46] stringr_1.5.0       splines_4.3.1       parallel_4.3.1     
[49] base64enc_0.1-3     vctrs_0.6.3         sandwich_3.1-0     
[52] SparseM_1.81        jsonlite_1.8.7      Formula_1.2-5      
[55] htmlTable_2.4.1     listenv_0.9.0       dr_3.0.10          
[58] foreach_1.5.2       gam_1.22-2          glue_1.6.2         
[61] parallelly_1.36.0   codetools_0.2-19    stringi_1.7.12     
[64] gtable_0.3.4        munsell_0.5.0       tibble_3.2.1       
[67] pillar_1.9.0        quantreg_5.95       htmltools_0.5.5    
[70] lava_1.7.3          R6_2.5.1            doParallel_1.0.17  
[73] evaluate_0.21       lattice_0.21-8      png_0.1-8          
[76] backports_1.4.1     MatrixModels_0.5-3  Rcpp_1.0.10        
[79] nlme_3.1-162        coda_0.19-4         gridExtra_2.3      
[82] checkmate_2.2.0     xfun_0.39           zoo_1.8-12         
[85] pkgconfig_2.0.3    

3. Folder Structure
-----------------------------------------------------------------------------------------------------
This folder contains the following files that can be used to reproduce all analysis and figures of the manuscript.

It contains two subfolders containing the following files:

./simulation/
	Main_simulation.R
	The main R script that performs the simulations and sources BasicFunctions.R, Table1.R, Table2.R,Table2(2).R, Table3.R,
	post_CRSIS.R, post_CRCSIS.R, FNR.R. Simulations was performed in parallel on several devices, mainly under window,
	but works with minor changes  under others. 

	Main_plot.R
	The main script to plot Figure 1 and Figure 2.

        BasicFunctions.R
        R code including loading required packages and  all the  functions defined for convenience .

	Table1.R
        R code including the function performing Example 5.1 simulations represented in Table 1.  

	Table2.R
        R code including the function performing Example 5.2 simulations represented in Table 2.

	Table2(2).R
        R code including the function performing Example 5.3 simulations represented in Table 2.   
  
	Table3.R
        R code including the function performing Example 5.4 simulations represented in Table 3.   

        post_CRSIS.R
        post_CRCSIS.R
        R code including the function performing variable selection results for Example 5.4 simulations represented in Table 4.  
           
        FNR.R
        R code including the function performing Example 5.5 simulations represented in Figure 2.  

./application:

	Main_data_anlaysis.R
	The main R script that performs the data analysis and sources BasicFunctions.R.  

        data.csv
	The public bladder cancer dataset of Dyrskjot et al. (2007) described in the manuscript (Section 6). The dataset contains the
	survival times of 404 patients and the gene expression measurements of 1381 genes for each patient, as well as clinical
	predictors such as age, sex, stage, WHO grade and treatment methos. 

        Table5.R
        R code  performing screening results in Table S1 and variable selection estimation results or CR-SIS, CR-CSIS, CR-SJS and
	crSIRS  in Table 5.

        Table6.R
        R code  performing C-index and IBS results in Table 6 for different penalized methos for  CR-SIS and CR-CSIS.

        Figure3.R
        R code generating prediction error curve estimates for three penalization methods and CR-SIS, CR-CSIS, CR-SJS and crSIRS
	and generating all subplots of Figure 3.
	  
 
4. Instruction on reproducing the results
-----------------------------------------------------------------------------------------------------

Step 1: Unzip all files into a folder and set the folder in previous step as the working directory in Rstudio. 

Step 2: To reproduce  simulation results presented in the manuscript, source files in subfolder ./simulation.
		 Main_simulation.R and Main_plots.R will output  results for Table 1 to 4, and Figure 1-2. 

        Note : The running time for CR-CSIS and CR-SJS might be a little longer, especially for CR-CSIS.

Step 3: To reproduce results for real data analysis, source files in subfolder ./application.
	Main_data_analysis.R will output results for Table 5, 6, S1 and Figure 3.
        
