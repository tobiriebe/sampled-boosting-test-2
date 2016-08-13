
library(BatchExperiments)
library(dplyr)

#
reg = loadRegistry("mytest-files")
showStatus(reg)

#getJobInfo(reg, findExpired(reg)) 
#getJobInfo(reg, findErrors(reg)) 
#getJobInfo(reg, findDone(reg)) 
