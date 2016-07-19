# Outputs for pMCMC time series.
pmcmcSummary <- function( pmcmcdata )
{
  par(mfrow = c(3, 1))
  plot(ts(pmcmcdata))
  acf(pmcmcdata)
  hist(ts(pmcmcdata), 30)

}

pmcmcSummary2d <- function( pmcmcdata )
{
  par(mfrow = c(4, 2))
  
  plot(ts(pmcmcdata[1:1]))
  plot(ts(pmcmcdata[2:2]))

  plot(rollmean(pmcmc_timeseries[1:1],10))
  plot(rollmean(pmcmc_timeseries[2:2],10))
  
  acf(pmcmcdata[1:1])
  acf(pmcmcdata[2:2])
  
  hist(ts(pmcmcdata[1:1]), 30)
  hist(ts(pmcmcdata[2:2]), 30)
  
}
setwd("../output")

pmcmcdata_timeseries = read.table("pmcmc_timeseries.txt")
pmcmcSummary2d(pmcmcdata_timeseries)

