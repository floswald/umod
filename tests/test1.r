# testing util_module

library(testthat)
library(umod)
setwd("~/git/umod/tests")
# setup objects

n = 5    # number of states
k = 5    # number of savings choices by state
m = 3    # number of discrete labor choices by state

cash   <- matrix(1:n,n,m)
cash   <- cash + matrix(0:2,n,m,byrow=TRUE)
labo   <- seq(from=0,to=1,length=m)
saving <- matrix(seq(from=0,to=8,length=k),n,k,byrow=TRUE)
EV     <- log(outer(1:n,1:n))
hsize  <- sample(0:2,size=n,replace=TRUE)
pars   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6,tau=1,myNA=-1000)
pars2   <- list(theta=0.2,phival=0.9,mu=0.6,gamma=1.4,cutoff=0.1,alpha=-0.6,tau=0.5,myNA=-1000)


res <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars)
res2 <- util_module(cashR=cash, saveR=saving, EVR=EV, hsizeR=hsize, laborR=labo, par=pars2)

context("testing that C output is sane")
test_that("res is a list",{
  expect_that( is.list(res), is_true() )})
test_that("res$values is a numeric matrix (n,m)",{
  expect_that( is.matrix(res$values) & all.equal(dim(res$values),c(n,m)) & !any(is.na(res$values)), is_true() )})
test_that("res$cons is a numeric matrix (n,m)",{
  expect_that( is.matrix(res$cons) & all.equal(dim(res$cons),c(n,m)) & !any(is.na(res$cons)), is_true() )})
test_that("res$saving is a numeric matrix (n,m)",{
  expect_that( is.matrix(res$saving) & all.equal(dim(res$saving),c(n,m)) & !any(is.na(res$saving)), is_true() )})
test_that("res$dchoiceL is matrix (n,1) of integers",{
  expect_that( is.matrix(res$dchoiceL) & all.equal(dim(res$dchoiceL),c(n,1)) & !any(is.na(res$dchoiceL)), is_true() )})
test_that("res$maxL is matrix (n,1) of doubles",{
  expect_that( is.matrix(res$maxL) & all.equal(dim(res$maxL),c(n,1)) & !any(is.na(res$maxL)), is_true() )})





# with NA values
# R function can't handle NAs

# for (i in 1:m){
#   tmp <- cash[,i]-saving2
#   
#   util <- ufun_labouR_disc(e=tmp,s=hsize,l=labo[i],params=pars)
#   W <- util + EV
#   rres[,i] <- apply(W,1,max)
#   rsave[,i] <- saving[cbind(1:n,apply(W,1,which.max))]
#   rcons[,i] <- tmp[cbind(1:n,apply(W,1,which.max))]
# }
# dchoiceL <- matrix(apply(rres,1,which.max),n,1)
#rsave <- matrix(rsave[cbind(1:n,dchoiceL)],n,1)
#rcons <- matrix(rcons[cbind(1:n,dchoiceL)],n,1)


# Rres <- list(values=rres,saving=rsave,cons=rcons,dchoiceL=dchoiceL,maxL=matrix(apply(rres,1,max),n,1))

# context("testing that C output is equal to R with NAs")
# test_that("R and C output are identical",{
#   expect_that( all.equal(Rres, res2),is_true()) })
# 


