
#test file for MT & CT
Rcpp::sourceCpp("/Users/orayasrim/Documents/MassTest/MassTesting/ShinyApp/ViralLoad.cpp")
source("/Users/orayasrim/Documents/MassTest/MassTesting/ShinyApp/buildFigures.R")
#plots (daily testing) another will have testing every 4 days for example ( noraml test period and outbreak test period)
#default 1 day and additional testing with 0 or 2 days ( tracing delay)
#fraction traced along x axis () ( probtraced given )
#prob detect symptoms (depends on disease) -> can set to 0 to ignore -> modify based on disease
#time to peak to symptoms -> modify based on disease
#days to peak -> depends on disease
#contact per day is 13, test delay = 10 hours, fraciso = 100 or 1, frac test = 0.9 to do thihnk about turn off contact tracing for people who dont test, 
#maxprobtranmit per exposure = 0.3, relative decline slope = 1, max days after peak. = 20 , mask effect = 0, minlogpcr = +2.5, initial log load = -2.5

# #base
# normalTestPeriod = 1 #adjust scenario
# outbreakTestPeriod = 1 #adjust scenario
# tracingDelay = 0 #adjust scenario
# fractionTraced = 0.1 #adjust this one
# probDetectSymptoms = 0 #adjust based on disease
# timeFromPeaktoSymptoms = 1 #adjust based on disease
# daysToPeak = 6 #adjust based on disease
# contactsPerDay = 13
# testDelay = 10
# fracIso = 1
# fracTest = 0.90
# maxProbTransmitPerExposure = 0.30
# relativeDeclineSlope = 1
# maxDaysAfterPeak = 20
# maskEffect = 0
# minLogPCRViralLoad = 2.5
# initialLogLoad = -2.5
# R0 = 3 # adjust based on disease 

#example diseases
# #influenza
#   probTest = .50 #adjust scenario
#   normalTestPeriod = 4 #adjust scenario
#   outbreakTestPeriod = 4 #adjust scenario
#   tracingDelay = 0 #adjust scenario
#   #fractionTraced = 0.1 #adjust this one
#   probDetectSymptoms = 0.693*probTest #adjust based on disease
#   timeFromPeaktoSymptoms = -2 #adjust based on disease
#   daysToPeak = 3 #adjust based on disease
#   contactsPerDay = 13
#   testDelay = 10
#   fracIso = 1
#   fracTest = 0.90
#   maxProbTransmitPerExposure = 0.30
#   relativeDeclineSlope = 1
#   maxDaysAfterPeak = 20
#   maskEffect = 0
#   minLogPCRViralLoad = 2.5
#   initialLogLoad = -2.5
#   R0 = 2.0 #adjust based on disease
  
# #SARS-CoV
# probTest = .80 #adjust scenario
# normalTestPeriod = 1 #adjust scenario
# outbreakTestPeriod = 1 #adjust scenario
# tracingDelay = 0 #adjust scenario
# fractionTraced = 0.1 #adjust this one
# probDetectSymptoms = 0.925*probTest #adjust based on disease
# timeFromPeaktoSymptoms = -5 #adjust based on disease
# daysToPeak = 7 #adjust based on disease
# contactsPerDay = 13
# testDelay = 10
# fracIso = 1
# fracTest = 0.90
# maxProbTransmitPerExposure = 0.30
# relativeDeclineSlope = 1
# maxDaysAfterPeak = 20
# maskEffect = 0
# minLogPCRViralLoad = 2.5
# initialLogLoad = -2.5
# R0 = 2.4 # adjust based on disease
  
 
test_period_list <- seq(0,20, by = 2)

disease_list <- c("1918 Influenza", "SARS-CoV-1")

fractionTraced_list <- seq(0, 0.9, by = 0.1)

dt_final = rbindlist(llply(test_period_list, function(test_period){
    
  normalTestPeriod = test_period
  outbreakTestPeriod = test_period
    
  rbindlist(llply(disease_list, function(disease_name){
    #normalTestPeriod = 6 #adjust scenario
    #outbreakTestPeriod = 6 #adjust scenario
    tracingDelay = 0 #adjust scenario
    #fractionTraced = 0.1 #adjust this one
    
    contactsPerDay = 13
    testDelay = 10
    fracIso = 1
    fracTest = 0.90
    maxProbTransmitPerExposure = 0.30
    relativeDeclineSlope = 1
    maxDaysAfterPeak = 20
    maskEffect = 0
    minLogPCRViralLoad = 2.5
    initialLogLoad = -2.5

    
    if(disease_name == "1918 Influenza"){
      #influenza
      probTest = .50 #adjust scenario & disease
      probDetectSymptoms = 0.693*probTest #adjust based on disease
      timeFromPeaktoSymptoms = -2 #adjust based on disease
      daysToPeak = 3 #adjust based on disease
      R0 = 2.0 #adjust based on disease
    }
    else{
      #SARS-CoV
      probTest = .80 #adjust scenario
      probDetectSymptoms = 0.925*probTest #adjust based on disease
      timeFromPeaktoSymptoms = -5 #adjust based on disease
      daysToPeak = 7 #adjust based on disease
      R0 = 2.4 # adjust based on disease
    }
    
    numOutbreaks = 120
    endDay = 120
    maxSize = 300
    
    rbindlist(llply(fractionTraced_list, function(fractionTraced){
      params = c(
        normalTestPeriod = normalTestPeriod*24  ,
        outbreakTestPeriod = outbreakTestPeriod*24 ,
        ContactTracingDelay = tracingDelay ,
        ProbTracedGivenInfectorDetected = fractionTraced,
        ProbDetectSymptoms = probDetectSymptoms,
        timeFromPeakToSymptoms = timeFromPeaktoSymptoms,
        timeToPeak = daysToPeak*24,
        maxTimeAfterPeak= 24*30,
        timeFromPeakTo0 = 24*5,
        contactsPerHour = contactsPerDay/24, testDelay = testDelay, fracIso = fracIso, fracTest = fracTest, 
        precision = 0.2, maxProbTransmitPerExposure = maxProbTransmitPerExposure,
        relativeDeclineSlope = relativeDeclineSlope, maxTimeAfterPeak = 24*maxDaysAfterPeak, 
        maskEffect = maskEffect, minLogPCRViralLoad = minLogPCRViralLoad, initialLogLoad = initialLogLoad, probTransmitMid = 8.9e6,logLimitOfDetection = 2.5, testPeriod = 1,fracQuar = 1 
        ,infectHParam = 0.51)
      
      
      #output final dt 
      params["logPeakLoad"] = computePeakViralLoad(params["timeToPeak"], targetR0 =R0, params)
      
      caseData = rbindlist(llply(1:numOutbreaks, function(i){
        dt = data.table(branchingModel(endDay = endDay, maxSize = maxSize, params)) 
        dt[,RunNumber := i]
        
        return(dt)
        }))
      re_filt <- caseData[hourNotInfectious<SimulationEndHour]
      re <- mean(re_filt$NumInfected)
      x = data.table(Re = re, Disease = disease_name, FractionTraced = fractionTraced, R0 = R0, testPeriod = test_period)
    return(x)
    }))
  }))
}))


p <- ggplot(data = dt_final[Disease == "1918 Influenza"]) + geom_line(aes(x = FractionTraced, y = Re),colour = "#00A2FF", size = 0.5)+  theme_minimal() + labs(y = "Effective Reproduction Number (Re)", x = "Proportion of Contacts Traced", title = "Re vs Proportion of Contacts Traced for 1918 Influenza Pandemic")

#ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_influenza.pdf", p, width = 7, height = 4,device = "pdf")

p2 <- ggplot(data = dt_final[Disease == "SARS-CoV-1"]) + geom_line(aes(x = FractionTraced, y = Re),colour = "#00A2FF", size = 0.5)+  theme_minimal() + labs(y = "Effective Reproduction Number (Re)", x = "Proportion of Contacts Traced", title = "Re vs Proportion of Contacts Traced for SARS-CoV-1")

#ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_SARS.pdf", p2, width = 7, height = 4,device = "pdf")


a <- ggplot(dt_final) +
  aes(
    x = FractionTraced,
    y = Re,
    colour = testPeriod,
    group = testPeriod
  ) +
  geom_line() +
  scale_color_distiller(palette = "Set2", direction = 1) +
  labs(
    x = "Proportion of Contacts Traced",
    y = "Effective Reproduction Number (Re)",
    title = "Re vs Proportion of Contacts Traced"
  ) +
  theme_minimal() +
  facet_wrap(vars(Disease))

# a <- ggplot(dt_final) +
#   aes(
#     x = FractionTraced,
#     y = Re,
#     colour = factor(testPeriod),
#     group = testPeriod,
#     fill = factor(testPeriod)
#   ) +
#   geom_line() + scale_colour_discrete(name = "testPriod") +
#   labs(
#     x = "x lab title ",
#     y = "y lab title",
#     title = "title main"
#   ) +
#   theme_minimal() +
#   facet_wrap(vars(Disease))


a <- ggplot(dt_final) +
  aes(
    x = FractionTraced,
    y = Re,
    colour = factor(testPeriod),
    group = testPeriod,
    fill = factor(testPeriod)
  ) +
  geom_line() + scale_colour_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(11), name = "Test Period (Days)") + 
  labs(
    x = "Proportion of Contacts Traced",
    y = "Effective Reproduction Number (Re)",
    title = "Re vs Proportion of Contacts Traced"
  ) +
  theme_minimal() + theme(legend.title=element_text(size=7),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size=7))+ 
  facet_wrap(vars(Disease)) 

ggsave("/Users/orayasrim/Documents/MassTest/Figures/re_vs_contact_trace_test_period.pdf", a, width = 7, height = 4,device = "pdf")

