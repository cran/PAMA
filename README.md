# Introduction of PAMA
'PAMA' implements Partition-Mallows model for rank aggregation where the rankers' quality are different. Rank aggregation aims to achieve a better ranking list given multiple observations. The problem of discerning reliability of rankers based only on the rank data is of great interest to many practitioners. By dividing the ranked entities into two disjoint groups, i.e., relevant and irrelevant/background ones, and incorporating the Mallows model for the relative ranking of relevant entities, PAMA can not only distinguish quality differences among the rankers but also provide the detailed ranking information for relevant entities.

This package provides both Bayesian inference and Maximum likelihood estimation (MLE). It can handle partial list as well. When covariates information is available, this package can make inference by incorporating the covariate information. More information can be found in the paper "Integrated Partition-Mallows Model and Its Inference for Rank Aggregation". The paper is accepted by Journal of the American Statistical Association.


# Structure of PAMA
Four main R functions for PAMA model.
1. PAMA.B.R  This function implements Bayesian inference of PAMA model.
2. PAMA.F.R  This function implements Maximum Likelihood estimation of PAMA model.
3. PAMA.PL.R This function implements Bayesian inference of PAMA model with partial lists.
4. PAMA.Cov.R This function implements Bayesian inference of PAMA model with covariates.

# Examples

```{r}
library(PAMA)
data<-NBANFL()
#To use function 'PAMA.B':
results<-PAMA.B(data$NBA,16,iter=10)
```



