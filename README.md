## Multivariate Galaxy module for PCA, PLS(-DA), and OPLS(-DA)

### Description

**Version:** 2.3.0  
**Date:** 2016-01-24  
**Author:** E. Thevenot (W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Licence:** CeCILL  

### Installation

* Configuration file: **multivariate_config.xml**
* Image files: 
    + **static/images/multivariate_workflowPositionImage.png**   
    + **static/images/multivariate_workingExampleImage.png**  
* Wrapper file: **multivariate_wrapper.R**  
* R packages  
    + **batch** from CRAN  
> install.packages("batch", dep=TRUE)  
    + **ropls** from Bioconductor  
> install.packages("batch", dep=TRUE)  
> source("http://www.bioconductor.org/biocLite.R")  
> biocLite("ropls")      

### News

##### CHANGES IN VERSION 2.3.0

NEW FEATURES  

    o **Predictions** now available (see the 'Samples to be tested' argument)  
    o OPLS(-DA): **Predictive and Orthogonal VIP** are now computed (see the 'comments' section)  
    o **Multiclass PLS-DA** implemented (see the 'comments' section)  
    
***