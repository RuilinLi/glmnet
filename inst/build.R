formatR::tidy_file('/Users/ruilinli/glmnet/R/allclasses.R')
tools::package_native_routine_registration_skeleton('/Users/ruilinli/glmnet', con='/Users/ruilinli/glmnet/src/glmnet_init.c')
devtools::document('/Users/ruilinli/glmnet')
install.packages('/Users/ruilinli/glmnet', repo=NULL, type='source')

library(glmnet)
a = list('/Users/ruilinli/plink-ng/toy_data.pgen')
