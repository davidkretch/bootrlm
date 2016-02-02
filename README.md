# bootrlm

bootrlm estimates robust linear models via bootstrap in R, implemented in C++ 
for speed. It currently supports MM-estimation, which is a robust estimation 
method that provides both high breakdown point and high efficiency. bootrlm 
is modeled after the rlm function in package MASS, some of whose code has been 
incorporated.

## Usage

```R
devtools::install_github("davidkretch/bootrlm")
library(bootrlm)
data(stackloss)
bootrlm_fit <- bootrlm(stack.loss ~ ., stackloss, r = 1000, method = "MM")
``` 