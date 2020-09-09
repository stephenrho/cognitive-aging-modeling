
# A Tutorial on Cognitive Modeling for Cognitive Aging Research

Here are materials accompanying the tutorial 'A Tutorial on Cognitive Modeling for Cognitive Aging Research' https://psyarxiv.com/qsnea/

- `models/` - contains the `stan` implementations of the signal detection theory models discussed in the tutorial
- `data/` - contains the data set used in the example analysis
- `simulate-data.R` - contains the code used to simulate the data set
- `fit_SDT.R` - uses the `rstan` package to fit the models
- `MPT` - contains the files needed to fit the multinomial processing tree model

To run the analysis you will need [R](https://www.r-project.org/) + [R Studio](https://www.rstudio.com/) and the following packages: `rstan, bridgesampling, loo, bayesplot, HDInterval`.

It can take a long time to fit the models. Therefore, to follow the tutorial without having to wait for `stan` to do the sampling, the fitted model objects can be downloaded from [here](https://drive.google.com/drive/folders/14gmtoYXKHMtZL7yjIzdmjKEhjIrmsGaq?usp=sharing). Once downloaded and unzipped, move the `.rds` files to the `models/` folder. They can then be loaded into `R` via:

```
name = readRDS("models/name.rds")
```
