# COLP

COLP is an R package for discovering causality from observational categorical data, developed and maintained by Yang Ni at Texas A&M University.

The package can also be downloaded at https://web.stat.tamu.edu/~yni/files/COLP_0.1.0.tar.gz.


#### Reference:  

Ni, Y. (2022). [Bivariate Causal Discovery for Categorical Data via Classification with Optimal Label Permutation](https://arxiv.org/pdf/2209.08579.pdf) *Advances in Neural Information Processing Systems (NeurIPS) 35*.


## A simple example
This example considers a real bivariate categorical data.
```{r dataload,echo=TRUE,  warning=FALSE, message=TRUE }
fit = COLP(CatPairs[[1]][[1]]$Diffwt,CatPairs[[1]][[1]]$Treat,algo="E")
fit$cd
```


