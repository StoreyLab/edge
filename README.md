dea: Differential Expression Analysis
====

Introduction
------
The edge package implements methods for carrying out differential 
expression analyses of genome-wide gene expression studies. Significance 
testing using the optimal discovery procedure and generalized likelihood 
ratio tests (equivalent to F-tests and t-tests) are implemented for general study 
designs. Special functions are available to facilitate the analysis of 
common study designs, including time course experiments. Other packages 
such as {\tt snm}, {\tt sva}, and {\tt qvalue} are integrated in edge to provide a wide range 
of tools for gene expression analysis.


### Installation and Documentation

To install, open R and type:

    install.packages("devtools")
    library("devtools")
    install_github("jdstorey/qvalue", build_vignettes = TRUE)
    install_github("jdstorey/dea", build_vignettes = TRUE)
    
Instructions on using edge can be viewed by typing:

    library("dea")
    browseVignettes("dea")

### Main functions
* build_models
* build_study
* odp
* lrt
* fit_models
* kl_clust
* apply_sva
* apply_snm
* apply_qvalue

### Quick start guide

To get started, first load the kidney dataset included in the package: 
```R
library(edge)
data(kidney)
names(kidney)
```
The kidney study is interested in determining differentially expressed genes with respect to age in kidney tissue. The `age` variable is the age of the subjects and the `sex` variable is whether the subjects were male or female. The expression values for the genes are contained in the `kidexpr` variable.
```R
kidexpr <- kidney$kidexpr
age <- kidney$age
sex <- kidney$sex
```

Once the data has been loaded, the user has two options to create the experimental models: `build_models` or `build_study`. If the experiment models are unknown to the user, `build_study` can be used to create the models:
```R
dea_obj <- build_study(data = kidexpr, adj.var = sex, tme = age, sampling = "timecourse")
full_model <- fullModel(dea_obj)
null_model <- nullModel(dea_obj)
```

The variable `sampling` describes the type of experiment performed, `adj.var` is the adjustment variable and `tme` is the time variable in the study. If the experiment is more complex then type `?build_study` for additional arguments.  

If the alternative and null models are known to the user then `build_models` can be used to make a deSet object:
```R
library(splines)
cov <- data.frame(sex = sex, age = age) 
null_model <- ~sex 
full_model <- ~sex + ns(age, df=4)
dea_obj <- build_models(data = kidexpr, cov = cov, null.model = null_model, full.model = full_model)
```

The `cov` is a data frame of covariates, the `null.model` is the null model and the `full.model` is the alternative model. The input `cov` is a data frame with the column names the same as the variables in the alternative and null models. Once the models have been generated, it is often useful to normalize the gene expression matrix using `apply_snm` and/or adjust for unmodelled variables using `apply_sva`.
```R
dea_norm <- apply_snm(dea_obj)
dea_sva <- apply_sva(dea_norm)

```

The `odp` or `lrt` function can be used on `dea_sva` to implement either the optimal discovery procedure or the likelihood ratio test, respectively:
```R
# optimal discovery procedure
dea_odp <- odp(dea_sva, bs.its = 30, verbose=FALSE)
# likelihood ratio test
dea_lrt <- lrt(dea_sva)
```

To access the proportional of null p-values estimate, p-values, q-values and local false discovery rates for each gene, use the function `qvalueObj`:
```R
qval_obj <- qvalueObj(dea_odp)
qvals <- qval_obj$qvalues
pvals <- qval_obj$pvalues
lfdr <- qval_obj$lfdr
pi0 <- qval_obj$pi0
```

See the vignette for more detailed explainations of the dea package.

