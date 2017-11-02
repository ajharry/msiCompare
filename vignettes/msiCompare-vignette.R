## ---- fig.show='hold',echo=FALSE-----------------------------------------
library(msiCompare)

## ---- fig.show='hold',echo=FALSE-----------------------------------------
library(Cardinal)

## ---- fig.show='hold'----------------------------------------------------
summary(s)

## ---- fig.show='hold'----------------------------------------------------
image(s$simSet, feature = 1)

## ---- fig.show='hold',cache=TRUE-----------------------------------------
fit_spautolm <- cass(msset = s$simSet, roiFactor = factor(s$simSet$diagnosis),
     logscale = F, thresholds = 1:5 )

fit_spautolm$results

