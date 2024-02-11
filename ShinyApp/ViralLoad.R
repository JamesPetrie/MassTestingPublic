require(data.table)
require(pracma)
require(plyr)

#require(pracma)

#folder = '/Users/orayasrim/Documents/MassTest/MassTesting/' #"/Users/orayasrim/Documents/MassTest/MassTesting/" # 
Rcpp::sourceCpp(paste0(folder, 'ViralLoad.cpp'))
# 
# 
# 
# # computes probability of infection for a typical contact given viral load
# probTransmit = function(viralLoad, params = NULL){
# 
#   return((viralLoad > 1)*(1 - exp(-0.2*(viralLoad^0.51)/(viralLoad^0.51 + (8.9e6)^0.51))))
#  }
# 
# 
# 
# 
# 
# # computes probability of a positive PCR result given viral load
# probPositive = function(viralLoad, params = NULL){
#   if(viralLoad > 1e3){return(1)}
#   else{return(0)}
# }
# 
# # compute viral load at a given time relative to infection
# # todo: decide if 10^0 is right intercept
# computeViralLoad = function(time, viralParams){
#   # if(time < 0){
#   #   logLoad = 0
#   # }else if(time <= viralParams$timeToPeak){
#   #   logLoad = viralParams$logPeakLoad*(time/viralParams$timeToPeak)
#   # }else if(time <= viralParams$timeToPeak + viralParams$timeFromPeakTo0){
#   #   logLoad = viralParams$logPeakLoad*(1 - ((time-viralParams$timeToPeak)/viralParams$timeFromPeakTo0))
#   # }else{
#   #   logLoad = 0
#   # }
#   
#   logVals = (time >= 0)*(time <= viralParams$timeToPeak)*viralParams$logPeakLoad*(time/viralParams$timeToPeak) +
#     (time > viralParams$timeToPeak)*(time <= viralParams$timeToPeak + viralParams$timeFromPeakTo0)*viralParams$logPeakLoad*(1 - ((time-viralParams$timeToPeak)/viralParams$timeFromPeakTo0))
#   return(10^logVals)
# }
# 
# 
# # computes expected transmissions over specified hour range relative to day of infection
# sumTransmissions = function(startTime, endTime, params){
#   val = integrate( lower = startTime, upper = endTime, f = function(tVals){
#     viralLoads = computeViralLoad(tVals, params)
#     return(params$contactsPerHour*probTransmit(viralLoads, params))
#   })
#   return(val$value)
# }
# 
# 
# #params = list(timeToPeak = 6*24, timeFromPeakTo0 = 8*24, logPeakLoad = 9, contactsPerHour = 13/24)
# 
# #sumTransmissions(0,24*12, params)
# 
# 
# # computes the fraction of transmissions that occur after receiving a positive test result
# fracAfterPositive = function(params){
#   # average over test timing offsets
#   
#   testPeriod = params$testPeriod
#   testDelay = params$testDelay
#   stopTime = params$timeToPeak + params$timeFromPeakTo0
#   
#   if(testPeriod <= 0) return(NA)
#   
#   # todo: fix so doesn't use same point twice (0 and testPeriod)
#   transmissionsAfter = mean(sapply(seq(0, testPeriod, length.out = 30), function(offset){
#     remainingProb = 1.0
#     expectedTransmissions = 0.0
#     testTime = offset
#     
#     # testTimes = seq(offset, offset + stopTime, by = testPeriod)
#     # probPos = sapply(testTimes, function(time){
#     #   probPositive(computeViralLoad(time, params), params)
#     # })
#     # 
#     # remainingProb = cumprod(1-probPos)
#     # 
#     # probEvent = shift(remainingProb, fill = 1)*probPos
#     # 
#     # expectedTransmissions = sapply(testTimes, function(time){
#     #   sumTransmissions(time + testDelay, stopTime, params)
#     # })
#     # 
#     # return(sum(probEvent*expectedTransmissions))
#     # 
#     # 
# 
#     while(testTime <= stopTime ){
#       probPos = probPositive(computeViralLoad(testTime, params), params)
# 
#       expectedTransmissions = expectedTransmissions + remainingProb*probPos*sumTransmissions(testTime + testDelay, stopTime, params)
# 
#       remainingProb = remainingProb*(1-probPos)
# 
#       testTime = testTime + testPeriod
#     }
#     return(expectedTransmissions)
#   })
#   
#   )
#   
#   return(transmissionsAfter/sumTransmissions(0, stopTime, params))
# }

# params = list(timeToPeak = 2*24, timeFromPeakTo0 = 6*24, logPeakLoad = 11, contactsPerHour = 13/24, testPeriod = 1*24, testDelay = 12)
# sumTransmissions(0,24*30, params)
# fracAfterPositive(params)
# 


# choose values for fracIso, fracTest, containment strategies
# iterate over timeToPeak 
# set timeFromPeakTo0 = 2*timeToPeak
# for each strategy:
  # use bisection to find logPeakLoad such that Re = 1
  # compute expected transmissions up to peak (1 day after peak?)
# plot (expected Transmissions before peak such that Re=1) vs timeToPeak for each strategy


evaluateStrategy = function(params){
  
  riseTimes = seq(36, 12*24, length.out = 5 + 30*params["precision"])
  dt = rbindlist(llply(riseTimes, function(timeToPeak){
    maxLoad = 16
    
    testParams = copy(params)
    testParams["logPeakLoad"] = maxLoad
    testParams["timeToPeak"] = timeToPeak
    testParams["timeFromPeakTo0"] = timeToPeak/testParams["relativeDeclineSlope"]
    

    maxR0 = sumTransmissionsStable(0,testParams["timeToPeak"] + testParams["timeFromPeakTo0"], testParams)
    if(maxR0 < 1){
      # disease controllable even without interventions - don't interpret these values
      return(data.table(MaxR0 = NA, Slope = testParams["logPeakLoad"]/timeToPeak, TimeToPeak = timeToPeak))
    }
    
    if(1 > sumTransmissionsStable(0,timeToPeak + testParams["timeFromPeakTo0"], testParams)*
       ((1-fracAfterPositive(testParams)*params["fracIso"]*params["fracTest"]))*
       (1- testParams["maskEffect"])){
      cutoffLoad = maxLoad
    }else{
      cutoffLoad = bisect(a = 0, b = maxLoad, maxiter = 10 + 30*params["precision"], fun = function(logPeakLoad){
        testParams["logPeakLoad"] = logPeakLoad
  
        
        R0 = sumTransmissionsStable(0,timeToPeak + testParams["timeFromPeakTo0"], testParams)
        if(R0 <= 0) return(1)
        
        fracAfter = fracAfterPositive(testParams)
        
        Re = R0*(1-fracAfter*params["fracIso"]*params["fracTest"])*(1- testParams["maskEffect"])
        return(1 - Re)
      })$root
    }
    
    testParams["logPeakLoad"] = cutoffLoad
    
    R0 = sumTransmissionsStable(0,testParams["timeToPeak"] + testParams["timeFromPeakTo0"], testParams)
    return(data.table(MaxR0 = R0, Slope = testParams["logPeakLoad"]/timeToPeak, TimeToPeak = timeToPeak))
  }))
  dt[ , TestDelay := params["testDelay"]]
  dt[, TestPeriod := params["testPeriod"]]
  dt[, FracIso := params["fracIso"]]
  dt[,FracTest := params["fracTest"]]

  return(dt)
}
# 
# inputParams = list(contactsPerHour = 13/24, testPeriods = 24*c(1,2,3,4), testDelay =  24, fracIso = 0.95, fracTest = 0.95)
# 
# 
# 
# 
# dt = rbindlist(llply(inputParams$testPeriods, function(testPeriod){
#   params = copy(inputParams)
#   params$testPeriod = testPeriod
#   evaluateStrategy(params)
# }))
# 
# 
# 
# ggplot(dt, aes(x = TimeToPeak/24, y = MaxR0, colour = paste("TestFreq=1/", TestPeriod/24))) + geom_line() + scale_y_continuous(breaks = 0:10) + scale_x_continuous(breaks = 0:10)

#p = ggplot()
#for(strat in input$strategies){
#  dt =
#}




