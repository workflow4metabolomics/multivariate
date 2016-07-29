## Multivariate analysis by PCA, PLS(-DA), and OPLS(-DA)
#### A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) project

### Description

**Version:** 2.3.2  
**Date:** 2016-05-15  
**Author:** Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)   
**Email:** [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Licence:** CeCILL
**Reference history:** [W4M00001a_sacurine-subset-statistics](http://galaxy.workflow4metabolomics.org/history/list_published), [W4M00001b_sacurine_complete](http://galaxy.workflow4metabolomics.org/history/list_published)     
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

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
> source("http://www.bioconductor.org/biocLite.R")  
> biocLite("ropls")      

### Tests

The code in the wrapper can be tested by running the **tests/multivariate_tests.R** in R  

### News

##### CHANGES IN VERSION 2.3.2

INTERNAL MODIFICATION  

    o Modification of the **multivariate_wrapper.R** file to handle **ropls** package versions of 1.3.15 and above (i.e. after switching to S4 classes)
    
***

##### CHANGES IN VERSION 2.3.0

NEW FEATURES  

    o **Predictions** now available (see the 'Samples to be tested' argument)  
    o OPLS(-DA): **Predictive and Orthogonal VIP** are now computed (see the 'comments' section)  
    o **Multiclass PLS-DA** implemented (see the 'comments' section)  
    
***