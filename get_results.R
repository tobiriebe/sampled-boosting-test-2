
library(BatchExperiments)
library(dplyr)

#
reg = loadRegistry("mytest-files")
showStatus(reg)

#
#getJobInfo(reg, findExpired(reg)) %>% select(time.started, time.done, time.running, memory) %>% View
#showLog(reg, 1)

# reduce results
res = reduceResultsExperiments(reg, ids = findDone(reg), fun = function(job, res) res$aggr)

# ggplot
g = res %>%
  ggplot() + geom_jitter(aes(x = algo, y = mmce.test.mean), colour = "red", size = 5, width = 0.3, alpha = 0.5) +
  facet_grid(~n.trees)  + ylab("MMCE") + ggtitle("GBM") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank())
g
