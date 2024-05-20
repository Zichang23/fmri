## Introduction

In this section, we will show you how to use our functions to calculate the number of edge for certain combination of frequencies and networks. 

### Data

The data we will use is `Caltech_0051475_rois_aal.1D`, which is a rs-fMRI dataset for a control subject from the Autism Brain Image Exchange (ABIDE) repository. For more information about this repository, please visit [ABIDE Preprocessed](preprocessed-connectomes-project.org/abide/).


### Functions

Here we will use two functions, `mconn2` and `mnet`.

`mconn2` is a function to calculate connectivity matrix across certain frequencies.

`mnet` is a function to calculate number of edges for all brain networks.

Below are all the arguments in `mconn2`. `x` is the dataset in format ".D". alpha is the significance level, the default value is 0.05. `s` is the number of stationary time series within the dataset `x`. `tt` is the time point for each stationary time series segment. `freq0` are the specific frequencies we would like to focus.

```{}
mconn2(x, alpha = 0.05, s, tt, freq0)
```

Below is the argument in `mnet`. `x` is the connectivity matrix, which can be get from `mconn2`. We will first apply `mconn2` to the dataset `Caltech_0051475_rois_aal.1D` to get a connectivity matrix for each specific  frequncy. Then we will plug the connectivity matrix into `mnet` to get the number of edge for each network.

```{}
mnet(x)
```

```{r}
#read in data
ctrl <- read.table("Caltech_0051475_rois_aal.1D", header = F)
head(ctrl[,1:8])
```
```{r}
#standardize the dataset ctrl (optional step)
library(dplyr)
ctrl1 <- ctrl %>% mutate_all(~(scale(.) %>% as.vector))
#select specific frequencies
freq0 <- seq(0.01,0.1,0.01)
freq0
```

```{r}
#apply functions to the data
library(astsa)
myconn <- mconn2(x=ctrl1, alpha=0.05, s=1, tt=41, freq0=freq0)
mynet <- mnet(myconn)
mynet
```
