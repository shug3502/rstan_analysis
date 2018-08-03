#just playing with different outputs of the model, predictions etc
library(rstan)
library(dplyr)
library(ggplot2)
source('make_fake_data.R')

data_nularge = make_fake_data(th=c(58,182),nu=0.95,sig=1.5)
data_nusmall = make_fake_data(th=c(58,364),nu=0.95,sig=1.5)
df1 <- data.frame(rna=data_nularge[20,],cellID=1:16,nu=rep('WT'))
df2 <- data.frame(rna=data_nusmall[20,],cellID=1:16,nu=rep('OE'))

g <- full_join(df1,df2) %>% ggplot(aes(cellID,rna,color=factor(nu))) + geom_line()
h <- full_join(df1,df2) %>% ggplot(aes(cellID,rna)) + facet_wrap(~nu) + geom_line()
print(g)
print(h)
