require(ggforce)
require(ggplot2)
require(cowplot)
require(scales)

require(data.table)
require(Rcpp)
require(plyr)
require(tidyr)
require(metR)
library(Cairo)
library(wesanderson)


folder = '~/MassTesting/ShinyApp/' #"/Users/orayasrim/Documents/MassTest/MassTesting/" # 
source(paste0(folder, "ViralLoad.R"))
#source("~/MassTesting/outbreakBranching.R")
#Rcpp::sourceCpp("ViralLoad.cpp")


typicalInitalLogLoad = -2.5
typicalPcrLogLod = 2.0
typicalAntigenLogLod = 5.4 # https://doi.org/10.1016/j.cmi.2022.08.006
infectiousMid = 8.9e6 # from Ke et al
typicalInfectH = 0.51 # from Ke et al
typicalPCRSlope = 6
typicalAntigenSlope = 1.3
typicalPCRMaxSens = 0.995
typicalAntigenMaxSens = 0.9 # https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2812586
#infectiousMid = 8.9e4 # test


defaultParams = c( maskEffect = 0,precision = 0.25, contactsPerHour = 13/24, maxProbTransmitPerExposure = 0.3,
                   relativeDeclineSlope = 1, maxTimeAfterPeak = 24*30, initialLogLoad = typicalInitalLogLoad, 
                   probTransmitMid = infectiousMid, infectHParam = typicalInfectH, timeFromPeakToSymptoms = 0, fracTransmitSymptoms = 1 ,
                   maxSensitivity = 0.995, testSlope = 6)


theme_set(theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) + theme(legend.background = element_rect(fill = "white")) + theme_half_open() + background_grid()  + 
            theme(text = element_text(size=20), axis.text = element_text(size=20)))





generateControllabilityFigure = function(testPeriods, inputParams){
  dt = rbindlist(llply(testPeriods, function(testPeriod){
    params = copy(inputParams)
    params["testPeriod"] = testPeriod
    evaluateStrategy(params)
  }))
  if(nrow(dt)>0){
    
    freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24), PeriodLabel = paste(TestPeriod/24)), by = TestPeriod]
    setkey(freqNames, by = "TestPeriod")
    freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
    freqNames[, PeriodLabel := factor(PeriodLabel, levels = PeriodLabel)]
    dt = merge(dt, freqNames, by = "TestPeriod")
    
    dt[,DaysToPeak := TimeToPeak/24]
    
    pathogenDt = data.table(Pathogen = c("SARS-CoV-2 (Wild-type)", "SARS-CoV-2 (Omicron)", "Influenza (1918)", "SARS-CoV-1", "Measles"),
                            DaysToPeak = c(5.0, 5.0, 3.5, 10, 10), R0 = c(2.5, 8, 2.0, 2.5, 15))
    pathogenDt = pathogenDt[, list(DaysToPeak = DaysToPeak + 0.5*c(1,0,-1,0), R0 = R0 + R0*0.2*c(0, 1,0,-1)), by = Pathogen ]
    
    p = ggplot() +
      geom_line(data= dt, aes(x = DaysToPeak, y = MaxR0, colour = PeriodLabel), linewidth = 1.4) + scale_y_log10(breaks = c(1,2,3,4,6,8,10,12,15, 20), limits = c(0.98,20)) + scale_x_continuous(breaks = 0:12) + 
      guides(colour=guide_legend(title="Days between tests")) + xlab("Days to Peak Viral Load") +ylab(expression("R"["0"] ))+
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, group = Pathogen, label = Pathogen), fill = "plum3",size = 0.00 ,label.fontsize = 14, show.legend = F,  lty = "blank")+  theme(legend.position="bottom") +
      labs(title = "Effect of Mass Testing",  subtitle = "Maximum controllable R0 for different testing strategies")
    #geom_ellipse(data = data.table(), aes(x0 = 5, y0 = 3, a = 1, b = 1, angle = 0), fill = "orange", alpha = 0.4)
    
    #geom_ribbon(aes(ymin = 0, ymax = MaxR0, x = TimeToPeak/24, fill = FreqLabel), alpha = 0.5)
  }
  return(p)
}


generate2TestControllabilityFigure = function(testPeriods, params){
  dt = data.table(expand.grid(TestPeriod = testPeriods, TestType = c("Antigen", "PCR")))
  
  dt[TestType == "Antigen", LOD := typicalAntigenLogLod]
  dt[TestType == "Antigen", TestDelay := 0]
  dt[TestType == "Antigen", TestSlope := typicalAntigenSlope]
  dt[TestType == "Antigen", MaxSensitivity := typicalAntigenMaxSens]
  dt[TestType == "PCR", LOD := typicalPcrLogLod]
  dt[TestType == "PCR", TestDelay :=8]
  dt[TestType == "PCR", TestSlope := typicalPCRSlope]
  dt[TestType == "PCR", MaxSensitivity := typicalPCRMaxSens]
  
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    newParams["logLimitOfDetection"] = dt[i,LOD]
    newParams["testSlope"] = dt[i,TestSlope]
    newParams["maxSensitivity"] = dt[i,MaxSensitivity]
    
    
    return(merge(dt[i,], evaluateStrategy(newParams)))
    
  }))
  
  
  if(nrow(dt)>0){
    
    freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24), PeriodLabel = paste(TestPeriod/24)), by = TestPeriod]
    setkey(freqNames, by = "TestPeriod")
    freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
    freqNames[, PeriodLabel := factor(PeriodLabel, levels = PeriodLabel)]
    dt = merge(dt, freqNames, by = "TestPeriod")
    
    dt[,DaysToPeak := TimeToPeak/24]
    dt[, TestType := factor(TestType, levels = c("PCR", "Antigen"))]
    
    pathogenDt = data.table(Pathogen = c("SARS-CoV-2 (Wild-type)", "SARS-CoV-2 (Omicron)", "Influenza (1918)", "SARS-CoV-1", "Measles"),
                            DaysToPeak = c(5.0, 5.0, 3.5, 10, 10), R0 = c(2.5, 8, 2.0, 2.5, 15))
    pathogenDt = pathogenDt[, list(DaysToPeak = DaysToPeak + 0.5*c(1,0,-1,0), R0 = R0 + R0*0.2*c(0, 1,0,-1)), by = Pathogen ]
    
    dt[TestPeriod/24 ==1 , TestStrategy := paste0(TestType, " Every Day", " (", TestDelay, " hour delay)")]
    
    dt[TestPeriod/24 != 1 , TestStrategy := paste0(TestType, " Every ", TestPeriod/24, " Days", " (", TestDelay, " hour delay)")]
    
    p = ggplot() +
      geom_line(data= dt, aes(x = DaysToPeak, y = MaxR0, colour = TestStrategy, linetype = TestStrategy), linewidth = 1.4) +
      scale_color_manual(values = c("red", "dodgerblue3", "red", "dodgerblue3")) + 
      scale_linetype_manual(values = c(2, 2, 1, 1)) +
      scale_y_log10(breaks = c(1,2,3,4,6,8,10,12,15, 20), limits = c(0.98,20)) + scale_x_continuous(breaks = 0:12) + 
      xlab("Time to Peak Viral Load [Days]") + ylab(expression("R"["0"] )) + guides(linetype=guide_legend(title="Test Strategy"), colour =guide_legend(title="Test Strategy")) + 
      theme(legend.key.width = unit(2,"cm"))+ 
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, group = Pathogen, label = Pathogen),fill = "plum3",size = 0.0 ,label.fontsize = 14, show.legend = F, lty = "blank")
    
    #guides(colour=guide_legend(title="Test period [Days]")) + guides(linetype=guide_legend(title="Test Type")) 
    
    #labs( subtitle = "Maximum controllable R0 for different testing strategies")
    #geom_ellipse(data = data.table(), aes(x0 = 5, y0 = 3, a = 1, b = 1, angle = 0), fill = "orange", alpha = 0.4)
    
    #geom_ribbon(aes(ymin = 0, ymax = MaxR0, x = TimeToPeak/24, fill = FreqLabel), alpha = 0.5)
  }
  return(p)
}




generateSymptomsControllabilityFigure = function(testPeriods, symptomTransmissionFractions , params){
  
  dt = data.table(expand.grid(TestPeriod = testPeriods, FracTransmitSymptoms = symptomTransmissionFractions))

  
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["fracTransmitSymptoms"] = dt[i,FracTransmitSymptoms]
    
    
    return(merge(dt[i,], evaluateStrategy(newParams)))
  }))
  
  
  if(nrow(dt)>0){
    
    freqNames = dt[, list(FreqLabel = paste("1 /", TestPeriod/24), PeriodLabel = paste(TestPeriod/24)), by = TestPeriod]
    setkey(freqNames, by = "TestPeriod")
    freqNames[, FreqLabel := factor(FreqLabel, levels = FreqLabel)]
    freqNames[, PeriodLabel := factor(PeriodLabel, levels = PeriodLabel)]
    dt = merge(dt, freqNames, by = "TestPeriod")
    
    dt[,DaysToPeak := TimeToPeak/24]
    
    pathogenDt = data.table(Pathogen = c("SARS-CoV-2 (Wild-type)", "SARS-CoV-2 (Omicron)", "Influenza (1918)", "SARS-CoV-1", "Measles"),
                            DaysToPeak = c(5.0, 5.0, 3.5, 10, 10), R0 = c(2.5, 8, 2.0, 2.5, 15))
    pathogenDt = pathogenDt[, list(DaysToPeak = DaysToPeak + 0.5*c(1,0,-1,0), R0 = R0 + R0*0.2*c(0, 1,0,-1)), by = Pathogen ]
    
    dt[TestPeriod/24 ==1 , TestStrategy := paste0("PCR ", " Every Day")]
    
    dt[TestPeriod/24 != 1 , TestStrategy := paste0("PCR ", " Every ", TestPeriod/24, " Days")]
    
    dt[, Symptoms := FALSE]
    dt[FracTransmitSymptoms != 1, Symptoms := TRUE]
    
    p = ggplot() +
      geom_line(data= dt, aes(x = DaysToPeak, y = MaxR0, colour = TestStrategy, linetype = Symptoms), linewidth = 1.4) +
      scale_color_manual(values = c("red", "dodgerblue3", "red", "dodgerblue3")) + 
      scale_y_log10(breaks = c(1,2,3,4,6,8,10,12,15, 20), limits = c(0.98,20)) + scale_x_continuous(breaks = 0:12) + 
      xlab("Time to Peak Viral Load [Days]") +ylab(expression("R"["0"] )) + 
      theme(legend.key.width = unit(2,"cm"))+ 
      geom_mark_ellipse(data = pathogenDt, aes(x= DaysToPeak, y = R0, group = Pathogen, label = Pathogen),fill = "plum3",size = 0.0 ,label.fontsize = 14, show.legend = F, lty = "blank") +
      guides(linetype=guide_legend(title="Symptoms"), colour =guide_legend(title="Test Strategy"))
   }
  return(p)
}


replaceParams = function(params, timeToPeak,  logPeakLoad){
  newParams = copy(params)
  newParams["logPeakLoad"] = logPeakLoad
  newParams["timeToPeak"] = timeToPeak
  newParams["timeFromPeakTo0"] = timeToPeak/params["relativeDeclineSlope"]
  return(newParams)
}

computePeakViralLoad = function(timeToPeak, targetR0,params){
  testParams = copy(params)
  testParams["timeToPeak"] = timeToPeak
  testParams["timeFromPeakTo0"] = timeToPeak/params["relativeDeclineSlope"]
  
  testParams["logPeakLoad"] = 20
  if(sumTransmissionsStable(0,timeToPeak + testParams["timeFromPeakTo0"], testParams) < targetR0){
    return(NA)
  }

  
  cutoffLoad = bisect(a = 0, b = 20, maxiter = 15, fun = function(logPeakLoad){
    testParams["logPeakLoad"] = logPeakLoad
    
    
    R0 = sumTransmissionsStable(0,timeToPeak + testParams["timeFromPeakTo0"], testParams)
    
    return(targetR0 - R0)
  })$root
  
  return(cutoffLoad)
  
}

wrangleFracAfterPositive = function(R0 , timeToPeak, testPeriod, lod, testDelay, params){
  newParams = copy(params)
  newParams["timeToPeak"] = timeToPeak
  newParams["timeFromPeakTo0"] = timeToPeak
  newParams["testPeriod"] = testPeriod
  newParams["logLimitOfDetection"] = lod
  newParams["testDelay"] = testDelay
  newParams["logPeakLoad"] = computePeakViralLoad(timeToPeak, targetR0 = R0 , params)
  
  return(fracAfterPositive(newParams))
}

plotFracReduction = function(params, R0 = 3, testPeriods = c(24, 72), timesToPeak = 24*c(3,6,9), n = 30, showPooled = FALSE){
  
  dt = data.table(expand.grid(TimeToPeak = timesToPeak, TestPeriod = testPeriods, LOD = seq(1,7, length.out = n), TestDelay = seq(0,48, length.out = n))) 
  
  dt[, FracAfterPositive :=  wrangleFracAfterPositive(R0, TimeToPeak, TestPeriod, LOD, TestDelay, params), by = list(TimeToPeak, TestPeriod, LOD, TestDelay) ]
  # for each of 2, 6, 10 peak time, peak VL of 8
  # for each of 1, 2, 3 day test period
  # for range of test delays
  # for range of LODs
  # for adherence = 1
  # expand.grid
  # compute fraction reduction for each row
  
  
  
  breaks = c(1,  0.99, 0.98, 0.96, 0.92, 0.84, 0.68, 0.36, 0)
  p = ggplot(dt, aes(TestDelay, LOD, z=  FracAfterPositive)) +scale_y_continuous(breaks = 1:6, labels = label_math(expr=10^.x) ) +
    geom_contour_filled(breaks = breaks, alpha = 0.6 ) + 
    
    #geom_point(aes(x = 0, y = 5), size = 3.5, colour = "blue")+
    
    #geom_hline(yintercept = 3, linetype = "dashed", colour = "blue", size = 1.0) + 
    #geom_textcontour(breaks = breaks, straight = T , position = "jitter"  ) + 
    geom_contour(aes(z = FracAfterPositive),  breaks = breaks, colour = "grey15", linetype = "dashed") + 
    geom_text_contour(aes(z = FracAfterPositive),  breaks = breaks,rotate = TRUE, nudge_x = 1.5, nudge_y = 0.1, skip = 0, colour = "black", label.placer = label_placer_fraction(frac = 0.99)) + 
    scale_fill_manual(values = terrain.colors(11)) + 
    #theme(legend.position = "none")  +
    guides(fill=guide_legend(title="Fraction Transmissions\nPrevented / Adherence"))+
    theme(legend.position = "bottom")+
    xlab("Test Delay (hours)") + ylab(expression("LOD"["50"] ~ "(copies/ml)"))  
  
  if(length(testPeriods) == 1){
    p = p + facet_wrap( ~paste("Days to Peak Viral Load:", TimeToPeak/24 ))+ theme(strip.background = element_blank())
  }else{
    p = p + facet_grid( paste("Test Period (days):",  TestPeriod/24)  ~ paste("Days to Peak Viral Load:", TimeToPeak/24 ))+ theme(strip.background = element_blank())
  }
  
  
  
  
  p = p +   annotate(geom="rect",  xmin = 2, xmax= max(dt$TestDelay), ymin = 1.5, ymax=3.3, fill="blue", alpha=0.2) + 
    annotate("text", x = max(dt$TestDelay) - 5, y = 3.2, label = "PCR")
  
  if(showPooled){
    p = p +   annotate(geom="rect",  xmin = 4, xmax= max(dt$TestDelay), ymin = 3.7, ymax=4.3, fill="blue", alpha=0.2) + 
      annotate("text", x = max(dt$TestDelay) - 14, y = 4.2, label = "10x Pooled PCR")
  }
  
  p = p +   annotate(geom="rect",  xmin = 0, xmax= 2, ymin = 5, ymax=7, fill="blue", alpha=0.2) + 
    annotate("text", x =8, y = 6, label = "Antigen\nLOD50")
 
  
  return(p)
}


generateReVsR0AndTime = function(params, lod = typicalPcrLogLod, testDelay = 8, n = 30, label = ""){
  dt = data.table(expand.grid(TimeToPeak = 24*seq(2,12, length.out = n), 
                              R0 = seq(0.5,12, length.out = n), #10^seq(0,log10(20), length.out = n), 
                              TestPeriod = 24, LOD = lod, TestDelay = testDelay)) 
  
  dt[, FracAfterPositive :=  wrangleFracAfterPositive(R0, TimeToPeak, TestPeriod, LOD, TestDelay, params), by = list(R0, TimeToPeak, TestPeriod, LOD, TestDelay) ]
  
  dt[, Re := R0*(1-FracAfterPositive)] 
  dt[, Label := label]
  return(dt)
}

plotReVsR0AndTime = function(dt){
  
  
 
  
  #breaks = c(1,  0.99, 0.98, 0.96, 0.92, 0.84, 0.68, 0.36, 0)
  breaks = seq(2,0.0, by = -0.1) #c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
  p = ggplot(dt, aes( x = TimeToPeak/24, y = R0, fill = Re)) +
    #scale_y_log10() +
    #geom_contour_filled(aes(z = Re), breaks = breaks, alpha = 0.6 ) + 
    
    geom_raster() + 
    geom_contour(aes(z = Re),  breaks = breaks, colour = "grey15", linetype = "dashed") + 
    #geom_text_contour(aes(z = FracAfterPositive),  breaks = breaks,  colour = "black")+ #, label.placer = label_placer_fraction(frac = 0.95)) + 
    #scale_fill_manual(values = terrain.colors(length(breaks), rev = TRUE)) +  
    scale_fill_gradient2(low = "green",high = "red" , mid = "yellow", midpoint = 1, name = expression("R"["e"]))+
    #scale_fill_gradient(low = "chartreuse", high = "darkorange")+
    #scale_fill_gradientn(colours = terrain.colors(10)) + 
    #scale_fill_gradientn(colours = terrain.colors(10)) + 
    xlab("Days to Peak Viral Load") + ylab(expression("R"["0"]))+  theme(strip.background = element_blank())  + 
    #guides(fill=guide_legend(title="Re"))+
    #theme(legend.position = "bottom")+
    facet_wrap(~Label)
  
  
  
  return(p)
}



plot3Trajectories = function(R0, timeToPeak, params){
  dt = data.table(expand.grid(TimeToPeak = timeToPeak, Time = 24*seq(0,16, length.out = 200))) 
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = R0, params), by = TimeToPeak]
  dt[ ,ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  dt[ ,DailyTransmissions :=24*params["contactsPerHour"]*symptomMult(Time, replaceParams(params, TimeToPeak, LogPeakLoad))*probTransmit(ViralLoad,params), by = list( Time, ViralLoad, LogPeakLoad, TimeToPeak) ]
  dt[ ,TestSensitivity :=probPositive(ViralLoad, params), by = list( ViralLoad, TimeToPeak) ]
  
  
  p1 = ggplot(dt, aes(x = Time/24, y = ViralLoad)) + geom_line(linewidth = 1.4)   + xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+ theme(text = element_text(size=12), axis.text = element_text(size=12))+
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    theme(strip.background = element_blank())  + ylab("Viral Load (copies / ml)")+ 
    labs(title="Viral Load Trajectory")
  
  
  p2 = ggplot(dt, aes(x = Time/24, y = TestSensitivity))  + geom_line(linewidth = 1.4) +  xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+theme(text = element_text(size=12), axis.text = element_text(size=12))+
    theme(strip.background = element_blank())  + ylab("Test Sensitivity")+ 
    labs(title="Test Sensitivity Trajectory")
  
  p3 = ggplot(dt, aes(x = Time/24, y = DailyTransmissions)) + geom_line(linewidth = 1.4) + 
    xlab("Day Since Infection" ) +theme(text = element_text(size=12), axis.text = element_text(size=12))+
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank(), strip.text.x = element_blank()) + ylab("Expected Daily Transmissions")+
    labs(title="Infectiousness Trajectory")
  
  return(list(p1,p2,p3))
}

plotTrajectories = function(params , R0){
  
  dt = data.table(expand.grid(TimeToPeak = 24*c(3,7,10), Time = 24*seq(0,16, length.out = 200))) 
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = R0, params), by = TimeToPeak]
  dt[ ,ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  dt[ ,DailyTransmissions :=24*params["contactsPerHour"]*probTransmit(ViralLoad,params), by = list( ViralLoad, TimeToPeak) ]
  dt[ ,TestSensitivity :=probPositive(ViralLoad, params), by = list( ViralLoad, TimeToPeak) ]
  
  
  peakNames = dt[, list(PeakLabel = paste("Days to Peak: ", TimeToPeak/24)), by = TimeToPeak]
  setkey(peakNames, by = "TimeToPeak")
  peakNames[, PeakLabel := factor(PeakLabel, levels = PeakLabel)]
  dt = merge(dt, peakNames, by = "TimeToPeak")
  
  
  dt[, LogViralLoad := log10(ViralLoad)]
  dtLong = melt(dt, id.vars = c("Time", "PeakLabel"), measure.vars = c("LogViralLoad", "TestSensitivity","DailyTransmissions"))
  #dtLong[, value:=value/max(value), by = variable]
  dtLong[variable == "LogViralLoad", variable := "Log Viral Load"]
  dtLong[variable == "TestSensitivity", variable := "Test Sensitivity"]
  dtLong[variable == "DailyTransmissions", variable := "Daily Transmissions"]
  
  #ggplot(dtLong, aes(x = Time/24, y = value)) + facet_wrap( ~   PeakLabel + variable , nrow = 3, dir = "v") + geom_line() + 
  #theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.y = element_blank()) + xlab("Day Since Infection" ) +
  # scale_x_continuous(breaks = seq(0,16, by = 2))+
  # theme(strip.background = element_blank()) + ggtitle("Example viral load trajectory, test sensitivity, expected transmissions")
  
  p1 = ggplot(dt, aes(x = Time/24, y = LogViralLoad)) + facet_wrap( ~   PeakLabel  , nrow = 1) + geom_line(linewidth = 1.4) + 
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x = element_blank()) +# xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank())  + ylab("Viral Load\n(log10 copies / ml)") 
  
  p2 = ggplot(dt, aes(x = Time/24, y = TestSensitivity)) + facet_wrap( ~   PeakLabel  , nrow = 1) + geom_line(linewidth = 1.4) + 
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.x = element_blank()) +# xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank(), strip.text.x = element_blank())  + ylab("Test\nSensitivity") + ylim(0,1) 
  
  p3 = ggplot(dt, aes(x = Time/24, y = DailyTransmissions)) + facet_wrap( ~   PeakLabel  , nrow = 1) + geom_line(linewidth = 1.4) + 
    xlab("Day Since Infection" ) +
    scale_x_continuous(breaks = seq(0,16, by = 2))+
    theme(strip.background = element_blank(), strip.text.x = element_blank()) + ylab("Expected Daily\nTransmissions")
  
  plot_grid(p1,p2,p3, ncol = 1,align = "v",rel_heights = c( 1,1,1.1)) 
}

plotInfectiousness = function(params){
  dt = data.table(ViralLoad = 10^seq(0, 12, length.out = 100))
  dt[ ,ExpectedTransmissions :=params["contactsPerHour"]*24*probTransmit(ViralLoad, params), by = ViralLoad ]
  
  p = ggplot(dt , aes(x = ViralLoad, y = ExpectedTransmissions)) + geom_line() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
    xlab("Viral Load (copies/ml)") + ylab("Expected Daily\nTransmissions") +
    theme(text = element_text(size=12), axis.text = element_text(size=12)) 
  #+ labs(title="Infectiousness vs. Viral Load")
  return(p)
}

plotTestSensitivity = function(pcrParams, antigenParams = NULL){
  dt = data.table(ViralLoad = 10^seq(0, 12, length.out = 100), TestType = "PCR")
  dt[ ,ProbPositive :=probPositive(ViralLoad, pcrParams), by = ViralLoad ]
  
  if(!is.null(antigenParams)){
    dt2 = data.table(ViralLoad = 10^seq(0, 12, length.out = 100), TestType = "Antigen")
    dt2[ ,ProbPositive :=probPositive(ViralLoad, antigenParams), by = ViralLoad ]
    dt = rbind(dt,dt2)
  }
 
  p = ggplot(dt , aes(x = ViralLoad, y = ProbPositive, group = TestType))  + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
    xlab("Viral Load (copies/ml)") + ylab("Test Sensitivity") + ylim(0, 1) +
    theme(text = element_text(size=12), axis.text = element_text(size=12)) + geom_line()
    #geom_textline(aes(
     #  label = TestType),hjust = 0.9)
  #+ labs(title="Test Sensitivity vs. Viral Load")
  return(p)
}

plotViralLoad = function(params){
  dt = data.table(Time = 24*seq(0, 12, length.out = 100))
  dt[ ,ViralLoad := computeViralLoad(Time, params), by = Time]
  
  
  p = ggplot(dt , aes(x = Time, y = ViralLoad)) + geom_line() + scale_y_log10(labels = c("V0"), breaks = 10^params["initialLogLoad"] ) + 
    scale_x_continuous(labels = c(expression("\u03c4"[p]), expression("\u03c4"[p] + "\u03c4"[r])), breaks = c(params["timeToPeak"], params["timeToPeak"] + params["timeFromPeakTo0"])) + 
    xlab("Time") + ylab("Log10 Viral Load") + 
    theme(text = element_text(size=12), axis.text = element_text(size=12)) 
  #+ labs(title="Test Sensitivity vs. Viral Load")
  return(p)
}


plotDailyTransmissionsWithoutTesting = function(params){
  dt = data.table(Time = 24*seq(0, 12, length.out = 100))
  dt[ ,ViralLoad := computeViralLoad(Time, params), by = Time]
  dt[ ,DailyTransmissions := 24*params["contactsPerHour"]*symptomMult(Time, params)*probTransmit(ViralLoad,params), by = Time ]
  ggplot(dt,aes(x = Time, y = DailyTransmissions)) + geom_line()  + xlab("Days Since Infection") + ylab("Daily Transmissions\nWithout Testing")
}



plotShapeSensitivity = function(params){
  # 4 panel plot: viral trajectories, default expected daily transmissions, 
  # fraction detected vs time, modified expected daily transmissions
  
  # 4 scenarios with R0 = 3: (2 days to peak, 6 days to peak) x (truncated after peak, symmetric)
  
  # iterate over scenarios, calculate viral load over time for each, transmission rate for every scenario and hour,
  # fraction discovered for every scenario and hour, modified transmissions for every scenario and hour
  
  dt = data.table(expand.grid(TimeToPeak = 24*c(4), Time = 24*seq(0,12, length.out = 200), Truncated = c(TRUE, FALSE))) 

  truncatedParams = copy(params)
  truncatedParams["relativeDeclineSlope"] = 100
  targetR0 = 3
  
  dt[Truncated == TRUE, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = targetR0, truncatedParams), by = TimeToPeak]
  dt[Truncated == FALSE, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = targetR0, params), by = TimeToPeak]
  
  
  dt[Truncated == TRUE, ViralLoad := computeViralLoad(Time, replaceParams(truncatedParams, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  dt[Truncated == FALSE, ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  p1 = ggplot(dt,aes(x = Time/24, y = log10(ViralLoad), group = paste(Truncated, TimeToPeak), linetype = Truncated)) + geom_line()   +xlab("Days Since Infection") + ylab("Viral Load (log10 copies/ml)")
  

  
  dt[ ,DailyTransmissions :=24*params["contactsPerHour"]*probTransmit(ViralLoad,params), by = list( ViralLoad, TimeToPeak, Truncated) ]
  
  p3 = ggplot(dt,aes(x = Time/24, y = DailyTransmissions, group = paste(Truncated, TimeToPeak), linetype = Truncated)) + geom_line()  +xlab("Days Since Infection") + ylab("Default Transmissions")
  
  dt[ Truncated == TRUE,FracDetect := fracDetected(Time,replaceParams(truncatedParams, TimeToPeak, LogPeakLoad)), by = list(TimeToPeak, LogPeakLoad) ]
  dt[ Truncated == FALSE,FracDetect := fracDetected(Time,replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(TimeToPeak, LogPeakLoad) ]
  
  
  p2 = ggplot(dt,aes(x = Time/24, y = FracDetect, group = paste(Truncated, TimeToPeak), linetype = Truncated)) + geom_line()  +xlab("Days Since Infection") + ylab("Probability Infection Detected")
  
  
  dt[, ModifiedTransmissions := (1-FracDetect)*DailyTransmissions]
  p4 = ggplot(dt,aes(x = Time/24, y = ModifiedTransmissions, group = paste(Truncated, TimeToPeak), linetype = Truncated)) + geom_line()  +xlab("Days Since Infection") + ylab("Modified Transmissions")
  
  print(dt[,sum(ModifiedTransmissions/24)/targetR0, by = paste(Truncated, TimeToPeak)])
  print(dt[,mean(LogPeakLoad/TimeToPeak), by = paste(Truncated, TimeToPeak)])
  
  p = plot_grid(p1 + theme(legend.position = "None"),p2 + theme(legend.position = "None"),p3 + theme(legend.position = "None") + ylim(0,3.5),p4 + theme(legend.position = c(0.6, 0.8))+ ylim(0,3.5), nrow = 2, labels = c("A","B","C","D"))
  
  return(p)
}


# todo: modify so that test frequency is varied, either PCR (LOD 3, Delay 12) or Antigen (LOD 5, delay 0) are used, and DaysToPeak in c(3,6,9)
# colour by days to peak, linetype by test type
# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotEffectTestFreq = function(R0, params){
  dt = data.table(expand.grid(TestType = c("Antigen", "PCR"), TestPeriod  = floor(24/seq(0.1,2, length.out = 24)), TimeToPeak = 24*c(3,6,9)))
  dt[TestType == "Antigen", LOD := typicalAntigenLogLod]
  dt[TestType == "Antigen", TestDelay := 0]
  dt[TestType == "PCR", LOD := typicalPcrLogLod]
  dt[TestType == "PCR", TestDelay :=8]
  
  #dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =4.5, params), by = TimeToPeak]
  #dt[, LogPeakLoad := 8]
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/params["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    newParams["logLimitOfDetection"] = dt[i,LOD]
    
    if(dt[i,TestType] == "Antigen"){
      newParams["maxSensitivity"] = typicalAntigenMaxSens
      newParams["testSlope"] = typicalAntigenSlope
    }else{
      newParams["maxSensitivity"] = typicalPCRMaxSens
      newParams["testSlope"] = typicalPCRSlope
    }
    
    newParams["logPeakLoad"] =  computePeakViralLoad(newParams["timeToPeak"], targetR0 = R0 , params) #dt[i, LogPeakLoad]
    
    fracAfter = fracAfterPositive(newParams)
    
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter)))
    
  }))
  
  
  dt[TestType == "Antigen", TestType := "Antigen\n(0 hour delay)"]
  
  dt[TestType == "PCR", TestType := "PCR\n(8 hour delay)"]
  
  
  dt[, TestType := factor(TestType, levels = c("PCR\n(8 hour delay)", "Antigen\n(0 hour delay)"))]
  
  
  p = ggplot(dt,aes(x = 24/TestPeriod, y = FracAfterPositive, colour = as.factor(TimeToPeak/24), linetype = TestType)) + geom_line(size = 0.9) +
    xlab("Tests Per Day") + ylab("Fraction Transmissions\nPrevented / Adherence") + 
    guides(colour=guide_legend(title="Days to Peak\nViral Load"), linetype = guide_legend(title = "Test Type")) + 
    scale_x_log10(breaks = c(1/8, 1/4, 0.5, 1,2), labels= c("1/8", "1/4", "1/2", "1", "2")) + theme(text = element_text(size=14), axis.text = element_text(size=14))
  
  
  #+scale_y_continuous(trans=logit_trans(), breaks = c(0.16, 0.5, 0.84, 0.92, 0.96, 0.98, 0.99, 0.995))
  
  return(p)
}


# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotEffectTestDelay = function(params, R0){
  dt = data.table(expand.grid(TestDelay = seq(0, 72, length.out = 24), TestPeriod  = 24,TestType = c("PCR", "Antigen"), TimeToPeak = 24*c(3,6,9)))
  # todo: add other updated antigen params
  
  #dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =4.5, params), by = TimeToPeak]
  #dt[, LogPeakLoad := 8]
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/params["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    
    if(dt[i,TestType] == "Antigen"){
      
      newParams["logLimitOfDetection"] = typicalAntigenLogLod
      newParams["maxSensitivity"] = typicalAntigenMaxSens
      newParams["testSlope"] = typicalAntigenSlope
    }else{
      newParams["logLimitOfDetection"] = typicalPcrLogLod
      newParams["maxSensitivity"] = typicalPCRMaxSens
      newParams["testSlope"] = typicalPCRSlope
    }
    
    
    newParams["logPeakLoad"] =  computePeakViralLoad(newParams["timeToPeak"], targetR0 = R0 , params) #dt[i, LogPeakLoad]
    
    fracAfter = fracAfterPositive(newParams)
    
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter)))
    
  }))
  

  dt[, TestType := factor(TestType, levels = c("PCR", "Antigen"))]
  
  p = ggplot(dt, aes(x = TestDelay, y = FracAfterPositive, colour =  as.factor(TimeToPeak/24), linetype = TestType)) + geom_line(size = 0.9)+ 
    guides(colour=guide_legend(title="Days to Peak\nViral Load"), linetype=guide_legend(title="Test Type\n(Daily Testing)")) + ylab("Fraction Transmissions\nPrevented / Adherence") +
    xlab("Test Delay (hours)") + theme(text = element_text(size=14), axis.text = element_text(size=14))
  
  return(p)
}

# todo: change to daily cost per person (in dollars)
# p = plotPrevalenceCost(c(1,3,7), c(variableTestCost = 2, isolationCost = 5000, fixedAnnualizedDailyTestCost = 1))
# ggsave(paste0(folder,"figures/prevalenceCost.pdf"), p, width = 7, height = 7,device = "pdf"
plotPrevalenceCost = function(testPeriods, params){
  # need cost per test
  # need cost per isolation
  # need test frequencies
  
  variableTestCost = params["variableTestCost"]
  isoCost = params["isolationCost"]
  fixedAnnualizedDailyTestCost = params["fixedAnnualizedDailyTestCost"]
  
  dt = data.table(expand.grid(DailyFracInfected = 10^seq(-7, -1, length.out = 200), TestPeriod  = testPeriods))
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    dailyInfections = dt[i, DailyFracInfected]
    testPeriod = dt[i,TestPeriod]
    
    
    
    testFreq = 1/testPeriod
    
    
    testCost = testFreq*variableTestCost
    
    
    gdpPerCapita = 70e3 # gdp per person in USA
    
    dailyTestCostPerPerson = testCost*365/gdpPerCapita
    dailyIsoCostPerPerson = dailyInfections*isoCost*365/gdpPerCapita
    
    dailyTotalCost = dailyTestCostPerPerson + dailyIsoCostPerPerson
    
    
    
    return(data.table(TotalCost = dailyTotalCost, TestCost = dailyTestCostPerPerson, IsoCost = dailyIsoCostPerPerson,  TestPeriod = testPeriod, FractionInfectedDaily = dailyInfections))
  } ))
  
  dt[,AnnualizedFixedCostPerPerson := fixedAnnualizedDailyTestCost/TestPeriod]
  
  dt[,PeriodDescription := paste0(TestPeriod, " ($",round(AnnualizedFixedCostPerPerson, digits = 3), " per person per year fixed cost)" )]
  
  dtLong = melt(dt, id.vars = c(
    "FractionInfectedDaily", "PeriodDescription"), measure.vars = c("TotalCost",  "IsoCost") )
  
  dtLong[ variable == "TotalCost", variable := "Total Cost"]
  dtLong[ variable == "IsoCost", variable := "Isolation Cost"]
  
  
  
  
  ggplot() +geom_line(data = dtLong, aes(x= FractionInfectedDaily, y = value, linetype = variable, group = paste(variable,PeriodDescription)),  linewidth = 2) + geom_line(data = dtLong[variable == "Total Cost"], aes(x= FractionInfectedDaily, y = value, colour = PeriodDescription), linewidth = 2) +
    scale_x_log10() +scale_y_log10(limits = c(0.002, 0.5), n.breaks = 8) + labs(x = "Daily Fraction of Population Infected", y = "Fraction of GDP") + 
    guides(colour = guide_legend(nrow = 2)) + theme(legend.position = c(0.1, 0.8)) +guides(colour=guide_legend(title="Test period [Days]"),linetype=guide_legend(title="Cost Type")) +
    ggtitle("Daily Testing and Isolation Cost During Outbreak")
  
}

# plotImportCost = function(){
#   
#   sampleNumGenerations = function(Re, overdispersion=0.1, initialSize = 1){
#     if(Re >= 1.0) return(NA)
#     numInfected = initialSize
#     
#     count = 0
#     totalInfected = 0
#     while(count < 1e4 && numInfected > 0){
#       totalInfected = totalInfected + numInfected
#       count = count + 1
#       numInfected = sum(rnbinom(numInfected, size = overdispersion, mu = Re ))
#     }
#     return(data.table(NumGenerations = count, TotalInfected = totalInfected ))
#   }
#   
#   
#   estimateOutbreakDuration = function(Re, initialSize, tau = 5, safetyMargin = 10){
#     numGen = mean(sapply(1:10000, function(x) {dt = sampleNumGenerations(Re = Re, initialSize = initialSize); return(dt$NumGenerations)}))
#     return(tau*numGen + safetyMargin)
#   }
#   #mean(sapply(1:10000, function(x) {dt = sampleNumGenerations(Re = 0.9, initialSize = 10); return(dt$TotalInfected)}))
#   
#   
#   steadyImportCost = function(importRate, popSize, Re, controlTime, dailyCost, infectionCost, detectionSize,  fracPopTarget, probExportPerInfected){
#     
#     numDailyImports = importRate*popSize
#     
#     if(is.na(controlTime)){
#       return(dailyCost + infectionCost*importRate*1/(1-Re))
#     }else{
#       outbreakSize = detectionSize/(1-Re)
#       outbreakCost = infectionCost*outbreakSize/popSize + fracPopTarget*controlTime*dailyCost
#       
#       if(probExportPerInfected*outbreakSize >= 1){
#         
#         return(dailyCost + infectionCost*importRate*1/(1-Re))
#       }
#       
#       
#       expectedNewOutbreaks = 1/(1 - probExportPerInfected*outbreakSize)
#       
#       
#       
#       
#       return(numDailyImports*outbreakCost*expectedNewOutbreaks) # todo: consider recursive outbreak seeding
#     }
#   }
#   
#   popSize = 60e6
#   fracDt = fread("~/MassTesting/FracTransmissions.csv")
#   R0 = 2
#   incentiveCost = 5
#   logisticCost = 8
#   poolSize = 24
#   pcrCost = 48
#   isoCost = 5e3
#   testFraction = 0.95
#   isoFraction = 0.9
#   testPeriod = 2
#   testDelay = 1
#   
#   dt = rbindlist(llply(10^seq(-11, -2, length.out = 200), function(importRate){
#     
#     transFraction = fracDt[TestPeriod == testPeriod & TestDelay == testDelay, FracAfterPositive]
#     
#     testFreq = 1/testPeriod
#     
#     Re = R0*(1 - testFraction*transFraction*isoFraction)
#     
#     if(Re > 1){
#       cost = NA
#       dailyTestCostPerPerson = NA
#       dailyIsoCostPerPerson = NA
#     }
#     else{
#       newInfectFrac = importRate*1/(1-Re) # fraction new infections (import and transmit per day) relative to total pop
#       
#       if(poolSize == 1){
#         testCost = testFreq*testFraction*(incentiveCost + logisticCost + pcrCost)
#       }else{
#         testCost = testFreq*testFraction*(incentiveCost + logisticCost + 1/poolSize*pcrCost + poolSize*newInfectFrac*testPeriod*pcrCost)
#       }
#       
#       gdpPerCapita = 70e3 # gdp per person in USA
#       
#       dailyTestCostPerPerson = testCost*365/gdpPerCapita
#       dailyIsoCostPerPerson = newInfectFrac*testFraction*isoCost*365/gdpPerCapita
#       
#       cost = dailyTestCostPerPerson + dailyIsoCostPerPerson
#     }
#     
#     
#     return(data.table(Cost = cost, TestCost = dailyTestCostPerPerson, IsoCost = dailyIsoCostPerPerson, Re = Re, TestPeriod = testPeriod, ImportRate = importRate, TestDelay = testDelay, FractionInfectedDaily = importRate/(1-Re)))
#     
#   } ))
#   
#   dt[, Strategy := "Continuous Testing"]
#   
#   
#   
#   
#   dtTimeTarget = copy(dt)
#   dtTimeTarget[, Strategy := "Temporal Targeting"]
#   dtTimeTarget[, FracPopTarget := 1]
#   dtTimeTarget[, ProbExportPerInfected := 0]
#   
#   dtSpatioTempTarget = copy(dt)
#   dtSpatioTempTarget[, Strategy := "Spatio-temporal Targeting (10k person divisions)" ]
#   dtSpatioTempTarget[, FracPopTarget := 10e3/60e6]
#   dtSpatioTempTarget[, ProbExportPerInfected := 0.02]
#   
#   
#   dtTarget = rbind(dtSpatioTempTarget, dtTimeTarget)
#   
#   
#   detectionSize = 10
#   dtTarget[ , LocalOutbreakDuration := estimateOutbreakDuration(Re, detectionSize, 5, 10), by = Re]
#   
#   dtTarget[, Cost := steadyImportCost(ImportRate, popSize, Re, LocalOutbreakDuration, TestCost, 5e3/70e3, detectionSize, FracPopTarget, ProbExportPerInfected), by = list(Re, LocalOutbreakDuration, TestCost, ImportRate, TestPeriod, TestDelay,FracPopTarget, ProbExportPerInfected )]
#   
#   dt = rbind(dt, dtTarget, fill = TRUE)
#   
#   dt[,TestPeriod := as.factor(TestPeriod)]
#   ggplot(dt[ TestPeriod ==2 & TestDelay == 1],aes(x= ImportRate, y = Cost, colour =Strategy)) + geom_line(linewidth = 2) + 
#     geom_vline(xintercept = 1e-5, linetype = "dashed", linewidth = 2) + annotate("text", x=6.5*1e-5, y=0.01, label="1000 cases \nper day (UK)", size = 6 ) +
#     geom_vline(xintercept = 1e-9, linetype= "dashed", linewidth = 2) +annotate("text", x=6.5*1e-9, y=0.01, label="0.1 cases \nper day (UK)", size = 6 )+
#     theme(legend.position = "bottom") +  guides(colour = guide_legend(nrow = 3)) +
#     scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + scale_y_log10(limits = c(0.001, 0.5), labels = trans_format("log10", math_format(10^.x)))+
#     labs(x = "Daily Imported Cases Relative to Population Size", y = "Daily Fraction of GDP") + 
#     theme(text = element_text(size=24), axis.text = element_text(size=24))
#     #+ylim(0,0.3)
#   
#   
#   
# }


plotOutbreaks = function(numOutbreaks = 10, endDay = 90, maxSize = 1000, params){
  caseData = rbindlist(llply(1:numOutbreaks, function(i){
    dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
    dt[,RunNumber := i]
    return(dt)
  }))
  
  
  dt = caseData[,list(DailyInfected = .N), by = list(Day = floor(InfectedHour/24),RunNumber ) ]
  dt = data.table(dt %>% complete(nesting(RunNumber), Day = seq(0, endDay, 1), fill = list(DailyInfected = 0)))
  setkeyv(dt, c("Day", "RunNumber"))
  dt[ , CumulativeInfected := cumsum(DailyInfected), by = RunNumber]
  
  
  
  
  #print(mean(caseData[DetectedHour < 1000*24,DetectedHour - InfectedHour]))
  #ggplot(caseData[DetectedHour < 1000*24], aes(x = DetectedHour - InfectedHour)) + geom_histogram()
  
  #ggplot(dt) + geom_line(aes(x = Day, y = CumulativeInfected, group = as.factor(RunNumber)), alpha = 0.4)  + geom_smooth(aes(x = Day, y = CumulativeInfected)) #+ scale_y_log10()
  ggplot(dt) + geom_line(aes(x = Day, y = DailyInfected, group = as.factor(RunNumber)), alpha = 0.4)  + geom_smooth(aes(x = Day, y = DailyInfected))
  
  # days of undetected infection per outbreak
  sumData = caseData[, sum(pmin(15*24, DetectedHour - InfectedHour, TracedHour - InfectedHour ))/24, by = RunNumber]; 
  print(mean(sumData$V1))
  ggplot(sumData, aes(x= V1)) + geom_histogram()
  
  # todo: days of elevated testing per outbreak
  
  # todo: cost of isolation and quarantine per outbreak
  
  # todo: number of infections per outbreak
}

plotMultiplePreventedTransmissions = function(params){
  newParams = copy(params)
  newParams["testDelay"] = 24
  newParams["testPeriod"] = 24*4
  p1 = plotPreventedTransmissions(newParams) + theme(text = element_text(size=14)) #+ ggtitle("Test Every 2 Days with 24 Hour Delay") 
  
  newParams["testDelay"] = 12
  newParams["testPeriod"] = 24*2
  p2 = plotPreventedTransmissions(newParams)  + theme(text = element_text(size=14)) # + ggtitle("Test Every Day with 8 Hour Delay")
  
  p = plot_grid(p1, p2, nrow= 1, labels = c("F", "G"))
  
}

# plot fraction transmissions after positive vs test delay for selected testing frequencies for each of the 3 viral load trajectories
plotFracTransmissionsAfterPositive = function(testPeriods,params, R0){
  dt = data.table(expand.grid(TestDelay = seq(0, 72, length.out = 12), TestPeriod  = testPeriods, TimeToPeak = 24*c(3,7,10)))
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 =R0, params), by = TimeToPeak]
  
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["logPeakLoad"] = dt[i, LogPeakLoad]
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/params["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    
    fracAfter = fracAfterPositive(newParams)
    
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter)))
    
  }))
  peakNames = dt[, list(PeakLabel = paste("Days to Peak: ", TimeToPeak/24)), by = TimeToPeak]
  setkey(peakNames, by = "TimeToPeak")
  peakNames[, PeakLabel := factor(PeakLabel, levels = PeakLabel)]
  dt = merge(dt, peakNames, by = "TimeToPeak")
  p = ggplot(dt, aes(x = TestDelay, y = FracAfterPositive, colour =  as.factor(TestPeriod/24))) + geom_line(linewidth = 1.4)  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.01)) + guides(colour=guide_legend(title="Test Period [Days]")) +
    xlab("Test Delay [Hours]") + ylab("Fraction Transmissions \n After Positive Test")+ 
    theme(legend.position = c(0.7, 0.2)) + theme(text = element_text(size=16)) + facet_wrap(~PeakLabel, nrow = 1)+
    theme(strip.background = element_blank())+  theme(legend.position="bottom") 
  return(p)
}

plotPreventedTransmissions = function(params){
  
  dt = data.table(TimeToPeak = params["timeToPeak"])
  dt[, LogPeakLoad := computePeakViralLoad(TimeToPeak, targetR0 = 4.5, params), by = TimeToPeak]
  dt = rbindlist(llply(1:nrow(dt), function(i){
    fracDiscovered = fracDiscoveredByHour( replaceParams(params, dt[i, TimeToPeak], dt[i,LogPeakLoad])) 
    hours = 0:(length(fracDiscovered) -1)
    return(cbind(data.table(Time = hours, FracDiscovered = fracDiscovered), dt[i,]))
  }))
  
  
  dt[ ,ViralLoad := computeViralLoad(Time, replaceParams(params, TimeToPeak, LogPeakLoad)), by = list(Time, TimeToPeak, LogPeakLoad)]
  
  dt[ ,DailyTransmissions := 24*params["contactsPerHour"]*symptomMult(Time, params)*probTransmit(ViralLoad,params), by = list( Time, ViralLoad, TimeToPeak) ]
  
  fracAdhere = params["fracTest"]*params["fracIso"]
  
  dt[, NonAdhereTransmissions := DailyTransmissions*(1-fracAdhere)]
  dt[, RemainingTransmissions := NonAdhereTransmissions + DailyTransmissions*(fracAdhere*(1-FracDiscovered))]
  
  dt[, Day:= Time/24]
  
  
  p =  ggplot(dt) + #geom_ribbon(aes(x = Time , ymin = RemainingTransmissions, ymax = HourlyTransmissions), fill = "grey", alpha = 0.6) + 
    geom_line(aes(x = Day , y = DailyTransmissions ),  linetype = "solid", colour = "black", linewidth = 1) +
    geom_ribbon(fill = "orange" , colour = "orange", alpha = 0.6 ,aes(x = Day ,ymin = 0,  ymax =  NonAdhereTransmissions)) + 
    geom_ribbon(fill = "purple" , colour = "purple", alpha = 0.6, aes(x = Day , ymin =  NonAdhereTransmissions, ymax = RemainingTransmissions)) + 
    xlab("Days Since Infection") + ylab("Expected Transmissions per Day") + geom_line(aes(x= Day, y= 1-(FracDiscovered)*params["fracTest"]), linetype= "dashed")
  
  return(p)
}


# randomly choose delay, lod50, test frequency, time to peak, R0
# compute frac prevented vs {P - ((vi -vt)/m - d)}
checkPattern = function(numPoints, params){
  dt = data.table(TestDelay = 72*runif(numPoints),
                  TestPeriod = 24*7*runif(numPoints),
                  LOD50 = 7.5*runif(numPoints) - 2.5,
                  TimeToPeak = 24*10*runif(numPoints) + 24,
                  R0 = 10*runif(numPoints) + 1,
                  InfectiousMid = 10^runif(numPoints, min = 5, max = 9) )

  
  dt = rbindlist(llply(1:nrow(dt), function(i){
    
    newParams = copy(params)
    newParams["logLimitOfDetection"] = dt[i, LOD50]
    newParams["timeToPeak"] = dt[i,TimeToPeak]
    newParams["timeFromPeakTo0"] = newParams["timeToPeak"]/newParams["relativeDeclineSlope"]
    newParams["testPeriod"] = dt[i,TestPeriod]
    newParams["testDelay"] = dt[i,TestDelay]
    newParams["infectiousMid"] = dt[i,InfectiousMid]
    
    newParams["logPeakLoad"] = computePeakViralLoad(newParams["timeToPeak"], targetR0 = dt[i, R0], newParams) 
    
    if(is.na(newParams["logPeakLoad"] )){
      return(NULL)
    }
    
    fracAfter = fracAfterPositive(newParams)
    
    slope =( newParams["logPeakLoad"] -  newParams["initialLogLoad"] )/newParams["timeToPeak"]
    
    vi = log10(newParams["infectiousMid"] ) - 1 # rough estimate
    vt = newParams["logLimitOfDetection"] + 1
    
    val = (vi - vt)/slope - newParams["testPeriod"] - newParams["testDelay"]
    
    fracTest = runif(10, min = 0.5, max = 1)
    isoEffect = runif(10, min = 0.5, max = 1)
    adherence = fracTest*isoEffect
    
    return(cbind(dt[i,], data.table(FracAfterPositive = fracAfter, ExtraTime = val, Adherence = adherence)))
    
  }))
  
  
  dt[ , Re := R0*(1-FracAfterPositive*Adherence)]
  
  ggplot(dt, aes(x = ExtraTime/24, y= (1-Adherence)*R0, colour = Re < 1)) + geom_point() + geom_hline(yintercept = 1) + geom_vline(xintercept = 0)
  
  
  
  
  #ggplot(dt, aes(x = ExtraTime/24, y = FracAfterPositive )) + geom_point(alpha = 0.4) + ylab("Fraction Transmissions\nPrevented") + 
   # xlab("(Time to Infectiousness) - (Time to detectable viral load)\n- (Test Period) - (Test Delay)")
    #xlab("(Early Detection Window) - (Test Period)\n- (Test Delay)")
  
 
  
}




generateCovidFracPrevented = function(params){
  covidParams = copy(params)
  covidParams["timeToPeak"] =5*24 # assuming 5 days to peak for Wild-type Covid-19
  covidParams["timeFromPeakTo0"] = covidParams["timeToPeak"]*covidParams["relativeDeclineSlope"]
  covidParams["logPeakLoad"] = computePeakViralLoad(covidParams["timeToPeak"], targetR0 = 2.5, covidParams) # assuming R0=2.5 for Wild-type Covid-19
  
  testScenarios = list(list(TestDelay = 0, TestPeriod = 24, LogLimitOfDetection  = typicalAntigenLogLod, Label = "Rapid Antigen"),
                       list(TestDelay = 8, TestPeriod = 24, LogLimitOfDetection  = typicalPcrLogLod, Label = "Fast PCR (8 Hours)"),
                       list(TestDelay = 72, TestPeriod = 24, LogLimitOfDetection = typicalPcrLogLod, Label = "Slow PCR (72 Hours)")
  )
  
  # iterate over 3 testing scenarios then iterate over frac adherence

  dt = rbindlist(llply(testScenarios, function(testScenario){
    newParams = copy(covidParams)
    newParams["logLimitOfDetection"] = testScenario$LogLimitOfDetection
    newParams["testDelay"] = testScenario$TestDelay
    newParams["testPeriod"] = testScenario$TestPeriod
     
    fracPrevented = fracAfterPositive(newParams)
    
    dt = rbindlist(llply(seq(0,1, by = 0.01), function(fracTest){
      dt = rbindlist(llply(c(0.5,0.95), function(isoEffect){
        fracAdherence = fracTest*isoEffect
        return(data.table(TestType = testScenario$Label, FracPrevented = fracPrevented*fracAdherence, 
                          FracCombinedAdherence = fracAdherence, IsoEffect = isoEffect, FracTest = fracTest))
       }))
    }))
    
    # dt = rbindlist(llply(seq(0,1, by = 0.01), function(fracAdherence){
    #   return(data.table(TestType = testScenario$Label, FracPrevented = fracPrevented*fracAdherence, FracCombinedAdherence = fracAdherence))
    # }))
    return(dt)
  }))

  dt[, IsoEffect := factor(IsoEffect, levels = c(0.95, 0.5))]
  # todo? make x axis fraction testing, linetype isolation effectiveness
  p1 = ggplot(dt,aes(x = FracTest, y = FracPrevented, colour = TestType)) +  geom_line(aes(linetype = IsoEffect))+
    xlab("Fraction Tested") + ylab("Fraction Transmissions Prevented") + theme( legend.position = c(0.1, 0.75)) +
    ylim(0,1) + #scale_y_continuous(breaks = seq(0,1, by = 0.2), labels = paste0(signif(100*seq(0,1, by = 0.2), 1), "%")) + 
    guides(colour=guide_legend(title="Test Type\n(With Daily Testing)"), linetype = guide_legend(title = "Isolation\nEffectiveness")) +
    theme(legend.title=element_text(size=14),legend.text=element_text(size=14)) 
  
  
  
  covidParams["logLimitOfDetection"] = typicalPcrLogLod
  covidParams["testDelay"] = 8
  # iterate over 3 adherence scenarios then iterate over test frequency (for fast PCR)
  

  dt = rbindlist(llply(floor(exp(seq(log(12), log(24*32), length.out = 80))), function(testPeriod){
    newParams = copy(covidParams)
    
    newParams["testPeriod"] = testPeriod
    
    fracPrevented = fracAfterPositive(newParams)
    
    dt = rbindlist(llply(c(0.1, 0.5, 0.9), function(fracTest){
      dt = rbindlist(llply(c(0.5, 0.95), function(isoEffect){
        fracAdherence = fracTest * isoEffect
        
        return(data.table(TestType = "PCR", TestPeriod = testPeriod, FracPrevented = fracPrevented*fracAdherence, FracCombinedAdherence = fracAdherence, IsoEffect = isoEffect, FracTest = fracTest))
      }))
    }))
    return(dt)
  }))
  
  dt[, IsoEffect := factor(IsoEffect, levels = c(0.95, 0.5))]
  
  p2 = ggplot(dt,aes(x = 24/TestPeriod, y = FracPrevented, colour = as.factor(FracTest) )) + geom_line(aes(linetype = IsoEffect))+
    xlab("Tests Per Day") + ylab("Fraction Transmissions Prevented") + ylim(0,1)+
    guides(colour=guide_legend(title="Fraction Tested \n(With Fast PCR)", reverse = TRUE),  linetype = guide_legend(title = "Isolation\nEffectiveness")) + 
    scale_x_log10(breaks = c(1/32, 1/16, 1/8, 1/4, 0.5, 1,2), labels= c("1/32", "1/16","1/8", "1/4", "1/2", "1", "2")) + 
    theme(legend.position = c(0.1, 0.75),  legend.title=element_text(size=14),legend.text=element_text(size=14)) 

  
  interventions = data.table(Intervention = c("Only schools and\nuniversities closed",  "Most nonessential\nbusinesses closed"), Effect = c(0.379,  0.266))
  
  
  p1 = p1 + geom_hline(data = interventions, aes(yintercept = Effect), linetype = "dashed", colour = "black") + geom_text(data = interventions, aes(x = 0,y = Effect,label = Intervention),  colour = "black", vjust = 0.5, hjust = 0, size = 3.5 )
  
  
  p2 = p2 + geom_hline(data = interventions, aes(yintercept = Effect), linetype = "dashed", colour = "black") + geom_text(data = interventions, aes(x = 1/32,y = Effect,label = Intervention),  colour = "black", vjust = 0.5, hjust = 0, size = 3.5 )
  
  
  my_colors <- RColorBrewer::brewer.pal(6, "Dark2")
  p = plot_grid(p1 + scale_colour_manual(values = my_colors[1:3]),p2+ scale_colour_manual(values = my_colors[4:6]), labels = c("A", "B"))
  return(p)
}

generateControllabilityScenarios = function(testPeriods, params){
  # Fig 4A
  p1 = generate2TestControllabilityFigure(testPeriods, combineParams(params , c( fracIso = 0.95, fracTest = 0.95))) 
  
  # Fig 4B
  p2 =  generate2TestControllabilityFigure(testPeriods, combineParams(params , c(  fracIso = 0.8, fracTest = 0.7)))
  
  
  p = plot_grid(p1+ theme(legend.position = "None") ,p2+ theme(legend.position = "None"), labels = c("A (90% Adherence)", "B (56% Adherence)"))
  grobs <- ggplotGrob(p1+guides(color = guide_legend(nrow = 2, title = "Test Strategy"), linetype = guide_legend(nrow = 2, title = "Test Strategy") ) +
                        theme(legend.direction = "horizontal",
                              legend.justification = "left",
                              legend.box.just = "bottom") )$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]] 
  p = plot_grid( legend, p , nrow = 2, rel_heights = c(0.12, 1))
  return(p)
}

combineParams = function(defaultParams, newParams){
  return(c(newParams, defaultParams[! names(defaultParams) %in% names(newParams)]))
}

generateReportFigures = function(folder){
  

 
  # Fig 1C 
  p = plotInfectiousness(defaultParams)
  ggsave(paste0(folder,"../figures/infectiousness.pdf"), p, width = 2.5, height = 2.5,device = "pdf")
  
  # Fig 1B
  p = plotTestSensitivity(combineParams(defaultParams , c(logLimitOfDetection = typicalPcrLogLod)), combineParams(defaultParams , c(logLimitOfDetection = typicalAntigenLogLod, maxSensitivity = typicalAntigenMaxSens, testSlope = typicalAntigenSlope)))
  ggsave(paste0(folder,"../figures/testSensitivity.pdf"), p, width = 2.5, height = 2.5,device = "pdf")
  
  
  # Figs 1A, 1D, 1E
  plots = plot3Trajectories(2.75, 5*24, combineParams( defaultParams, c(logLimitOfDetection = 2, timeFromPeakToSymptoms = 0, fracTransmitSymptoms = 0.5)))
  
  ggsave(paste0(folder,"../figures/viralLoad.pdf"), plots[[1]], width = 4, height = 4,device = "pdf")
  ggsave(paste0(folder,"../figures/testSensitivityVsTime.pdf"), plots[[2]], width = 4, height = 4,device = "pdf")
  ggsave(paste0(folder,"../figures/infectiousnessVsTime.pdf"), plots[[3]], width = 4, height = 4,device = "pdf")
  
  # Fig 1F, 1G
  p = plotMultiplePreventedTransmissions(combineParams(defaultParams , c( testDelay = 8, fracIso = 0.9, fracTest = 0.9, logLimitOfDetection = 2, timeToPeak = 6*24, relativeDeclineSlope = 1)))
  ggsave(paste0(folder,"../figures/preventedTransmissions.pdf"), p, width = 10, height = 6,device = "pdf")
  
  # Fig 2
  p = generateCovidFracPrevented(combineParams(defaultParams, c(precision = 0.6,  timeFromPeakToSymptoms = 0.0, fracTransmitSymptoms = 0.5, precision = 5.0 )))
  ggsave(paste0(folder,"../figures/covidPreventedTransmissions.pdf"), p, width = 10, height = 6,device = "pdf")
  
  # Fig 2 Appendix - more symptoms
  p = generateCovidFracPrevented(combineParams(defaultParams, c(precision = 0.6,  timeFromPeakToSymptoms = -24, fracTransmitSymptoms = 0.25 , precision = 2.0)))
  ggsave(paste0(folder,"../figures/covidPreventedTransmissionsMoreSymptoms.pdf"), p, width = 10, height = 6,device = "pdf")
  
  
  # Fig 3A
  p1 = plotFracReduction(R0  = 3 ,testPeriods = 24*c(1), timesToPeak = 24*c(3,6,9), n = 20, showPooled = FALSE,combineParams(defaultParams ,c())) 
    
  # Fig 3B
  p2 = plotEffectTestDelay(R0 = 3, params= combineParams(defaultParams , c( precision = 1.0)))
  
  # Fig 3C
  p3 = plotEffectTestFreq(R0 = 3, params= combineParams(defaultParams , c( precision = 1.0)))
  
  # Fig 3
  p = plot_grid(p1, plot_grid(p2,p3, nrow = 1, labels = c( "B", "C")), nrow = 2, labels = c("A", "") , rel_heights = c(1.4,1))
  ggsave(paste0(folder,"../figures/testParams.pdf"), p, width = 10, height = 9,device = "pdf")
  
  # Fig4
  p = generateControllabilityScenarios(24*c(1,3), combineParams(defaultParams , c( testDelay = 8, precision = 0.4)))
  ggsave(paste0(folder,"../figures/controllabilityComparison.pdf"), p, width = 14, height = 8,device = "pdf")
  
  # Fig4 appendix - more symptoms
  p = generateSymptomsControllabilityFigure(24*c(1,3), c(1.0, 0.25), combineParams(defaultParams , c( testDelay = 8, fracIso = 0.95, fracTest = 0.95, timeFromPeakToSymptoms  = 0,  logLimitOfDetection = typicalPcrLogLod, precision = 0.4))) 
  #p = generateControllabilityScenarios(24*c(1,3), combineParams(defaultParams , c( testDelay = 8, timeFromPeakToSymptoms = -24, fracTransmitSymptoms = 0.25)))
  ggsave(paste0(folder,"../figures/controllabilityComparisonMoreSymptoms.pdf"), p, width = 14, height = 8,device = "pdf")
  
  # Fig4 appendix - lower Km
  p = generateControllabilityScenarios(24*c(1,3), combineParams(defaultParams , c( testDelay = 8, probTransmitMid = infectiousMid/100, precision = 0.4)))
  ggsave(paste0(folder,"../figures/controllabilityComparisonLowKm.pdf"), p, width = 14, height = 8,device = "pdf")
  
  # appendix fig: dependence of Re on R0 and time to peak
  dt1 = generateReVsR0AndTime(combineParams(defaultParams , c(precision = 0.5)) , lod = typicalPcrLogLod, testDelay = 8, n = 30, label = "PCR (8 hour delay)")
  
  dt2 = generateReVsR0AndTime(combineParams(defaultParams , c( precision = 0.5,maxSensitivity = typicalAntigenMaxSens, testSlope = typicalAntigenSlope)) , lod = typicalAntigenLogLod, testDelay = 0, n = 30, label = "Antigen")
  dt = rbind(dt1,dt2)
  p = plotReVsR0AndTime(dt)
  ggsave(paste0(folder,"../figures/ReVsR0AndTimeToPeak.pdf"), p, width = 10, height = 6,device = "pdf")
  
 
 # m = as.matrix(cast(dt2, R0 ~  TimeToPeak, value = "Re" ))
  #fig <- plot_ly(x = sort(unique(dt2$TimeToPeak/24)), y = sort(unique(dt2$R0)), z = m, ) %>% add_surface() %>%layout(scene = list(xaxis = list(title = "Days to Peak"), yaxis = list(title = "R0"), zaxis = list(title = "Re")))
  #m = as.matrix(cast(dt1, R0 ~  TimeToPeak, value = "Re" ))
  #fig <- plot_ly(x = sort(unique(dt2$TimeToPeak/24)), y = sort(unique(dt2$R0)), z = m, ) %>% add_surface() %>%layout(scene = list(xaxis = list(title = "Days to Peak"), yaxis = list(title = "R0"), zaxis = list(title = "Re")))
  
  
  #p = plot_grid(p2 + theme(legend.position = "none"), p1, nrow = 2, labels = c( "Antigen","PCR (8 hour delay)") )
  # Fig 3 Appendix version: more test periods
  #p = plotFracReduction(testPeriods =24*c(1,2,4), timesToPeak = 24*c(3,6,9), n = 10, showPooled = TRUE, combineParams(defaultParams , c(logPeakLoad = 8)))
  #ggsave(paste0(folder,"figures/fracReduction2DSupplement.pdf"), p, width = 10, height = 12,device = "pdf")
  

  #p = generateSymptomsControllabilityFigure(24*c(1,3), c(1.0, 0.25), combineParams(defaultParams , c( testDelay = 8, fracIso = 0.95, fracTest = 0.95, timeFromPeakToSymptoms  = -24,  logLimitOfDetection = typicalPcrLogLod))) 
  # Fig 3A Appendix = lower Km
  p = plotFracReduction(R0 = 3, testPeriods = 24*c(1), timesToPeak = 24*c(3,6,9), n = 20, showPooled = FALSE,combineParams(defaultParams ,c(probTransmitMid = infectiousMid/100))) 
  ggsave(paste0(folder,"../figures/fracReduction2DLowerKm.pdf"), p, width = 10, height = 6,device = "pdf")
  
 
 p = plotShapeSensitivity( combineParams(defaultParams , c(timeFromPeakToSymptoms = 0, fracTransmitSymptoms = 1.0, testPeriod = 24, testDelay = 8, logLimitOfDetection = typicalPcrLogLod,precision = 2 )))
 ggsave(paste0(folder,"../figures/viralLoadShapeFastTest.pdf"), p, width = 10, height = 10,device = "pdf")
 
 p = plotShapeSensitivity( combineParams(defaultParams , c(timeFromPeakToSymptoms = 0, fracTransmitSymptoms = 1.0, testPeriod = 48, testDelay = 24, logLimitOfDetection = typicalPcrLogLod, precision = 2 )))
 ggsave(paste0(folder,"../figures/viralLoadShapeSlowTest.pdf"), p, width = 10, height = 10,device = "pdf")
 
}


computeValuesForMainText =function(){
  # compute fraction reduction in transmissions for PCR with 8 hour delay for 
  # {3.5 days to peak viral load, R0=2.5} and {5 days to peak viral load, R0=2.5}
  fracReducInfluenza = wrangleFracAfterPositive(R0 = 2.5 , timeToPeak = 3.5*24, testPeriod = 24, lod = typicalPcrLogLod, testDelay = 8, defaultParams)
  fracReducCovid = wrangleFracAfterPositive(R0 = 2.5 , timeToPeak = 5*24, testPeriod = 24, lod = typicalPcrLogLod, testDelay = 8, defaultParams)
  
  print(paste0("Fraction reduction in Influenza transmissions with perfect adherence, daily PCR testing with 8 hour delay: ",fracReducInfluenza  ))

  print(paste0("Fraction reduction in Covid transmissions with perfect adherence, daily PCR testing with 8 hour delay: ",fracReducCovid ))
  # compute adherence needed for PCR with 8 hour delay and {5 days to peak viral load, R0=2.5} 
  # to prevent 1.5/2.5 and 0.5/1.5 transmissions
  
  # 1 = 2.5 *  (1 - fracReduc *adherence)
  # adherence = (1- 1/(2.5 ))/fracReduc
  highAdherence = (1- 1/(2.5 ))/fracReducCovid
  print(paste("Adherence needed to reduce Re from 2.5 to 1 with daily PCR testing with 8 hour delay:", highAdherence ))
  
  lowAdherence = (1- 1/(1.5 ))/fracReducCovid
  print(paste("Adherence needed to reduce Re from 1.5 to 1 with daily PCR testing with 8 hour delay:", lowAdherence ))
  
  
  fracReducCovid2DayTest = wrangleFracAfterPositive(R0 = 2.5 , timeToPeak = 5*24, testPeriod = 48, lod = typicalPcrLogLod, testDelay = 8, defaultParams)
  fracReducCovid4DayTest = wrangleFracAfterPositive(R0 = 2.5 , timeToPeak = 5*24, testPeriod = 96, lod = typicalPcrLogLod, testDelay = 8, defaultParams)
  print(paste("relative fraction reduction for testing every 2 days vs 4 days:", fracReducCovid2DayTest/fracReducCovid4DayTest))
  
  
  # compute Re for PCR with 72 hour delay and {5 days to peak viral load, R0=2.5}, with two previous adherence values
  
  fracReducCovidSlowTest = wrangleFracAfterPositive(R0 = 2.5 , timeToPeak = 5*24, testPeriod = 24, lod = typicalPcrLogLod, testDelay = 72, defaultParams)
  
  print(paste("Re for Covid daily PCR testing with 72 hour delay with high adherence:", 2.5 * (1 - highAdherence* fracReducCovidSlowTest )))
  
  print(paste("Re for Covid daily PCR testing with 72 hour delay with low adherence:", 2.5 * (1 - lowAdherence* fracReducCovidSlowTest )))
}


#inputParams = c(contactsPerHour = 13/24, testDelay = 24, fracIso = 0.9, fracTest = 0.9, precision = 0.15, maxProbTransmitPerExposure = 0.3, relativeDeclineSlope = 1.0, maxTimeAfterPeak= 24*30, logPeakLoad = 10, initialLogLoad = -2, logLimitOfDetection = 3)

#generateControllabilityFigure(c(24, 48), inputParams)
#plotTrajectories(inputParams)
#plotImportCost()

#plotFracTransmissionsAfterPositive(c(24, 48), inputParams)


#plotPreventedTransmissions(c(inputParams, c(testPeriod = 48, timeToPeak = 100, timeFromPeakTo0 = 96)))

