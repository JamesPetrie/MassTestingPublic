#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(shinycssloaders)



folder = ""  #'/Users/jpetrie/MassTesting/' #"/Users/orayasrim/Documents/MassTest/MassTesting/" #
source(paste0(folder, "buildFigures.R"))




ui <- fluidPage(useShinyjs(),
  h2("Controlling a Respiratory Virus Pandemic with Mass Testing"),

                sidebarLayout(
                  sidebarPanel(
                    p(h3("Parameters"),
                      actionButton("resetAll", "Reset all parameters")),
                  
                    # todo: put reset button here
                    div(id = "params",tabsetPanel(
                    
                    tabPanel("Interventions",
                      withMathJax(p('$$R_e = R_0 \\cdot (1 - \\gamma  \\cdot \\beta \\cdot \\sigma) \\cdot (1-\\lambda)$$')),
                      withMathJax(p('\\(R_e\\) is the effective reproduction number (average number of transmissions per case). If \\(R_e < 1\\) then the disease is under control.')),
                      
                      #sliderInput("contactsPerDay","Contacts Per Day:", min = 0,max = 40,value = 13),
                      sliderInput("fracTest",withMathJax(p('\\(\\gamma\\): Fraction of (homogeneous) population testing regularly:')), min = 0,max = 0.995,value = 0.90),
                      sliderInput("fracIso",withMathJax(p('\\(\\beta\\):  Isolation effectiveness (fraction reduction in transmissions for detected positives):')), min = 0,max = 0.995,value = 0.90),
                      withMathJax(p('\\(\\sigma\\):   Fraction of counterfactual transmissions occurring after receiving a positive test result. The calculation is shown in the Model tab (depending on viral load trajectory, test frequency, and test delay)')),
                      sliderInput("testDelay", "Test Delay [hours]:",min = 0, max = 72,value = 12),
                      checkboxGroupInput("testPeriods", "Days Between Tests", choices = list(1, 2,3,5,7,10,30), selected= (list(1)), inline = TRUE),
                      sliderInput("maskEffect",withMathJax(p('\\(\\lambda\\):  Fraction of transmissions prevented by masks:')), min = 0,max = 0.995,value = 0.0),
                      
                      
                    ),
                    
                    tabPanel("Test Sensitivity",
                             plotOutput("TestSensitivity", height="140px"),
                             withMathJax(p('$$S_V(V) =  S_{max} \\frac{1}{1 + e^{-k(\\text{log}_{10}(V)-\\text{log}_{10}(\\text{LOD}_{50}))}}$$')),
                             sliderInput("logLimitOfDetection",withMathJax(p('$$LOD_{50}  \\text{  (log10 copies / ml):}$$')), min = 0.0,max = 8.0,value = 2.0, step = 0.25), 
                             sliderInput("maxSensitivity",withMathJax(p('$$S_{max}  \\text{ (peak test sensitivity:)}$$')), min = 0.0,max = 1.0,value = 0.995, step = 0.005),
                             sliderInput("testSlope",withMathJax(p('$$k  \\text{ (controls width of intermediate region):}$$')), min = 0.0,max = 10.0,value = 6.0, step = 0.5),
                    ),
                    tabPanel("Infectiousness",
                             plotOutput("Infectiousness", height="140px"),
                             withMathJax(p('$$T_V(V) = N_c \\cdot(1 - e^{-\\theta \\frac{V^h}{V^h + K_m^h}})$$')),
                             sliderInput("maxProbTransmitPerExposure",withMathJax(p('$$\\theta \\text{ (max probability transmission per exposure):}$$')), min = 0.1,max = 0.9,value = 0.3), # 0.3 would be consistent with 95% of measles household contacts infected -> 1 - 0.7^8 = 0.94
                             sliderInput("contactsPerDay",withMathJax(p('$$N_c \\text{ (contacts per day}):$$')), min = 1,max = 50,value = 13), 
                             sliderInput("infectHParam",withMathJax(p('$$h \\text{ (controls width of transition region):}$$')), min = 0,max = 1,value = 0.51),
                             shinyWidgets::sliderTextInput("probTransmitMid",withMathJax(p('$$k_m \\text{ (midpoint of infectiousness):}$$')),
                                                           choices=c("10^3", "10^3.5", "10^4", "10^4.5", "10^5","10^5.5","10^6","10^6.5","10^7","10^7.5","10^8","10^8.5","10^9"),
                                                           selected= "10^7", grid = TRUE) 
                    ),
                    tabPanel("Viral Load",
                             plotOutput("ViralLoad", height="140px"),
                         #     withMathJax(p('$$\\text{log}_{10}(V(t)) = 
                         # \\begin{cases}
                         # \\displaylines{
                         # \\text{log}_{10}(V_0) \\cdot (1 - \\frac{t}{\\tau_p}) + \\text{log}_{10}(V_{p}) \\cdot \\frac{t}{\\tau_p}, & \\text{for } 0 < t \\leq \\tau_p \\\\
                         # \\text{log}_{10}(V_p) \\cdot (1 - \\frac{t - \\tau_p}{\\tau_r}) + \\text{log}_{10}(V_0) \\cdot \\frac{t - \\tau_p}{\\tau_r}, & \\text{for }  \\tau_p > t > \\tau_p + \\tau_r \\\\
                         # -\\infty, & \\text{otherwise }}
                         # \\end{cases}$$')),
                         withMathJax(p('$$\\tiny{
                              \\text{log}_{10}(V(t)) = 
                                \\begin{cases}
                                \\displaylines{
                              \\text{log}_{10}(V_0) \\cdot (1 - \\frac{t}{\\tau_p}) + \\text{log}_{10}(V_{p}) \\cdot \\frac{t}{\\tau_p}, & \\text{for } 0 < t \\leq \\tau_p \\\\
                              \\text{log}_{10}(V_p) \\cdot (1 - \\frac{t - \\tau_p}{\\tau_r}) + \\text{log}_{10}(V_0) \\cdot \\frac{t - \\tau_p}{\\tau_r}, & \\text{for }  \\tau_p > t > \\tau_p + \\tau_r \\\\
                              -\\infty, & \\text{otherwise }}
                              \\end{cases}}$$')),
                         
                             sliderInput("initialLogLoad",withMathJax(p('$$\\text{log}_{10}(V_0) \\text{ (Initial viral load (log10 copies / ml))}:$$')), min = -4, max = 0,value = -2.5, step = 0.25),
                             sliderInput("relativeDeclineSlope",withMathJax(p('$$\\frac{\\tau_p}{\\tau_r}$$')), min = 0.05,max = 5.0,value = 1.0),
                             sliderInput("maxDaysAfterPeak",withMathJax(p('$$\\text{maximum days after peak \n viral load infection ends:}$$')), min = 0,max = 20,value = 20)
                             
                             # todo: computation of expected transmissions after postive test
                             
                             # sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
                             
                    ),
                    # tabPanel("Symptoms",
                    #          plotOutput("DailyTransmissions", height="140px"),
                    #          #     withMathJax(p('$$\\text{log}_{10}(V(t)) = 
                    #          # \\begin{cases}
                    #          # \\displaylines{
                    #          # \\text{log}_{10}(V_0) \\cdot (1 - \\frac{t}{\\tau_p}) + \\text{log}_{10}(V_{p}) \\cdot \\frac{t}{\\tau_p}, & \\text{for } 0 < t \\leq \\tau_p \\\\
                    #          # \\text{log}_{10}(V_p) \\cdot (1 - \\frac{t - \\tau_p}{\\tau_r}) + \\text{log}_{10}(V_0) \\cdot \\frac{t - \\tau_p}{\\tau_r}, & \\text{for }  \\tau_p > t > \\tau_p + \\tau_r \\\\
                    #          # -\\infty, & \\text{otherwise }}
                    #          # \\end{cases}$$')),
                    #          sliderInput("fracDetectSymptoms",withMathJax(p('Fraction that notice symptoms')), min = 0.0,max = 1.0,value = 0.5),
                    #          
                    #          sliderInput("daysFromPeakToSymptoms",withMathJax(p('Days from peak viral load to symptoms')), min = -6,max = 6,value = 0),
                    #          
                    #          sliderInput("fracTransmitSymptoms",withMathJax(p('Fraction of transmissions that still occur after symptoms ')), min = 0,max = 1,value = 0.5)
                    #        
                    #         
                    #          
                    #          # todo: computation of expected transmissions after postive test
                    #          
                    #          # sliderInput("simPrecision","Simulation Precision:", min = 0,max = 1,value = 0.2)
                    #          
                    # )
                    ))),                       
                  
                  
                  # Show a plot of the generated distribution
                  mainPanel(
                    
                    tabsetPanel(
                      tabPanel("Intervention Effectiveness",
                               #HTML("<b>Scenario</b>: Novel airborne pandemic with no available vaccines or treatments. Global elimination unlikely so expecting imported cases<br/>"),
                               #HTML("<b>Goal</b>: Reduce the number of infections while waiting for vaccines<br/>"),
                               #HTML("<b>Challenge</b>: Because of the high cost of existing Non-Pharmaceutical Interventions (NPIs), many countries may be unwilling or unable to control the epidemic<br/>"),
                               #HTML("<b>Proposed solution</b>: Frequent (inexpensive and convenient) saliva PCR testing for most of the population. Generous financial support for isolation of people who test positive.<br/><br/><br/>"),
                               shinycssloaders::withSpinner(plotOutput("controlRegion", height="550px"), hide.ui = FALSE),
                               withMathJax(p('The graph is created by finding the largest \\(R_0\\) (modified by varying peak viral load) for each testing strategy such that \\(R_e \\leq 1\\) (described in intervention parameters tab). Pathogens below each line can be independently controlled by that testing strategy.')),
                               
                      ),
                      
                      
                      tabPanel("Model Details",
                               # todo: figure of characteristic curve
                               
                               h3("Example Viral Load, Infectiousness, and Test Sensitivity Trajectories (with \\(R_0=3\\))"),
                               shinycssloaders::withSpinner(plotOutput("Trajectories", height="500px"), hide.ui = FALSE),
                               h3("Fraction of transmissions after a positive test"),
                               #Todo: Add figures for fraction of transmissions occuring after positive test for each trajectory
                               p("\\(\\sigma\\), the fraction of counterfactual transmissions occuring after receiving a positive test  can be computed for each testing strategy and viral load trajectory."),
                               shinycssloaders::withSpinner(plotOutput("FracAfterPositive", height="300px"), hide.ui = FALSE),
                               h3("Calculation"),
                               p("Let \\(E[T(x)|\\pi]\\) be the expected number of transmissions on day x after infection, conditional on parameters \\(\\pi\\) describing the viral load trajectory.
             Let \\(P(Test(x) | \\pi)\\) be the probability that a test taken on day \\(x\\) is positive, also conditional on the viral load trajectory. 
             Define the function ProbAllNegative(x, offset, period) as the probability that all samples collected on or before day x are negative, 
             with the time between sequential tests set by the period variable, and the timing relative to infection set by the offset variable."),
                               withMathJax(p('$$ProbAllNegative(x, offset, period) = \\prod_{k=0}^{(x-offset)/period}(1-P(Test(k*period + offset) | \\pi))$$')),
                               p('The expected fraction of transmissions after a positive test is computed by averaging over test timing offset, with offset ~ unif(0,period):'),
                               withMathJax(p('$$\\sigma(testDelay, testPeriod) = \\frac{E_{offset}[\\sum_{j=0}^{\\infty} E[T(j)|\\pi] \\cdot (1 - ProbAllNegative(j-testDelay, offset, testPeriod))]}{\\sum_{j=0}^{\\infty} E[T(j)|\\pi]}$$')),
                               
                               
                               
                               
                      )
                      
                    )
                  )
                  
                  
                  #Todo: show peak image with peak viral load and time to peak, describe how figure generated
                )
)



 
 


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  observeEvent(input$resetAll, {
    reset("params")
  })
  #output$exampleImplementation <- 
  

  
   output$controlRegion <- renderPlot({
 
     probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
     
     #input$contactsPerDay
     # inputParams = c(contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, probTransmitMid = probTransmitMidNumeric, infectHParam = input$infectHParam,
     #                 precision = 0.1, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
     #                 relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
     #                 maskEffect = input$maskEffect, logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, 
     #                 fracTransmitSymptoms = input$fracTransmitSymptoms, timeFromPeakToSymptoms = 0, maxSensitivity = 0.995, testSlope = 6)
     # 
     inputParams = c(contactsPerHour = input$contactsPerDay/24, testDelay = input$testDelay, fracIso = input$fracIso, fracTest = input$fracTest, probTransmitMid = probTransmitMidNumeric, infectHParam = input$infectHParam,
                     precision = 0.1, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure,
                     relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                     maskEffect = input$maskEffect, logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, 
                     fracTransmitSymptoms = 1, timeFromPeakToSymptoms = 0, maxSensitivity = input$maxSensitivity, testSlope = input$testSlope)
     
     generateControllabilityFigure(24*as.numeric(input$testPeriods), inputParams)
     
     
    })
   
   output$Trajectories <- renderPlot({
     #reactive({
       #req(getIncperMedianlogContour()) # can't plot it until these values have been calculated
     probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
     # inputParams = c( logPeakLoad = 10, contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = probTransmitMidNumeric,infectHParam = input$infectHParam,
     #                  relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
     #                  logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, precision = 0.15, 
     #                  fracTransmitSymptoms = input$fracTransmitSymptoms, timeFromPeakToSymptoms = 0,maxSensitivity = 0.995, testSlope = 6)
     # 
     inputParams = c( logPeakLoad = 10, contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = probTransmitMidNumeric,infectHParam = input$infectHParam,
                      relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                      logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, precision = 0.15, 
                      fracTransmitSymptoms = 1, timeFromPeakToSymptoms = 0,maxSensitivity = input$maxSensitivity, testSlope = input$testSlope)
     
      plotTrajectories(inputParams, R0 = 3) 
       
    # })
   })
   
   output$FracAfterPositive <- renderPlot({
     probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
     # inputParams = c(  contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = probTransmitMidNumeric,infectHParam = input$infectHParam,
     #                  relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
     #                  logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, precision = 0.15, 
     #                  fracTransmitSymptoms = input$fracTransmitSymptoms, timeFromPeakToSymptoms = 0, maxSensitivity = 0.995, testSlope = 6)
     # 
     inputParams = c(  contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = probTransmitMidNumeric,infectHParam = input$infectHParam,
                       relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, 
                       logLimitOfDetection = input$logLimitOfDetection, initialLogLoad = input$initialLogLoad, precision = 0.15, 
                       fracTransmitSymptoms = 1, timeFromPeakToSymptoms = 0, maxSensitivity = input$maxSensitivity, testSlope = input$testSlope)
     
     plotFracTransmissionsAfterPositive(24*as.numeric(input$testPeriods), R0 = 3, inputParams)
   })
   
   # output$Infectiousness <- renderPlot({
   #   inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, 
   #                   relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     
     output$Infectiousness <- renderPlot({
       probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
       inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, probTransmitMid = probTransmitMidNumeric, infectHParam = input$infectHParam,
                       relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak)
     plotInfectiousness(inputParams) 
   })
   # output$TestSensitivity <- renderPlot({
   #   probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
   #   inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, logLimitOfDetection = input$logLimitOfDetection, probTransmitMid = probTransmitMidNumeric,infectHParam = input$infectHParam,
   #                   maxSensitivity = 0.995, testSlope = 6)
     output$TestSensitivity <- renderPlot({
       probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
       inputParams = c(contactsPerHour = input$contactsPerDay/24, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure, relativeDeclineSlope = input$relativeDeclineSlope, maxTimeAfterPeak = 24*input$maxDaysAfterPeak, logLimitOfDetection = input$logLimitOfDetection, probTransmitMid = probTransmitMidNumeric,infectHParam = input$infectHParam,
                       maxSensitivity = input$maxSensitivity, testSlope = input$testSlope)
       
     plotTestSensitivity(inputParams) 
   })
   output$ViralLoad <- renderPlot({
     probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
     timeToPeak = 3.5*24
     timeFromPeakTo0 = timeToPeak/input$relativeDeclineSlope

     inputParams = c(logPeakLoad = 8 , maxTimeAfterPeak = 24*input$maxDaysAfterPeak, initialLogLoad = input$initialLogLoad,
                     timeToPeak =  timeToPeak, timeFromPeakTo0 = timeFromPeakTo0)
     plotViralLoad(inputParams) 
   })
   output$DailyTransmissions <- renderPlot({
     probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
     timeToPeak = 5*24
     timeFromPeakTo0 = timeToPeak*input$relativeDeclineSlope

     # inputParams = c(logPeakLoad = 8 , maxTimeAfterPeak = 24*input$maxDaysAfterPeak, initialLogLoad = input$initialLogLoad,
     #                 timeToPeak =  timeToPeak, timeFromPeakTo0 = timeFromPeakTo0, fracTransmitSymptoms = input$fracTransmitSymptoms, fracDetectSymptoms = input$fracDetectSymptoms, timeFromPeakToSymptoms = input$daysFromPeakToSymptoms*24,
     #                 contactsPerHour = input$contactsPerDay/24,  fracIso = 0, fracTest = 0, probTransmitMid = probTransmitMidNumeric, infectHParam = input$infectHParam,
     #                 precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure)
     inputParams = c(logPeakLoad = 8 , maxTimeAfterPeak = 24*input$maxDaysAfterPeak, initialLogLoad = input$initialLogLoad,
                     timeToPeak =  timeToPeak, timeFromPeakTo0 = timeFromPeakTo0, fracTransmitSymptoms = 1, fracDetectSymptoms = input$fracDetectSymptoms, timeFromPeakToSymptoms = 0,
                     contactsPerHour = input$contactsPerDay/24,  fracIso = 0, fracTest = 0, probTransmitMid = probTransmitMidNumeric, infectHParam = input$infectHParam,
                     precision = 0.2, maxProbTransmitPerExposure = input$maxProbTransmitPerExposure)
     
     plotDailyTransmissionsWithoutTesting(inputParams) 
   })

   # output$PrevalenceCost <- renderPlot({
   #   probTransmitMidNumeric = eval(parse(text = input$probTransmitMid))
   #   inputParams = c(variableTestCost = input$variableTestCost, isolationCost = input$isolationCost, fixedAnnualizedDailyTestCost = input$fixedAnnualizedDailyTestCost)
   #   plotPrevalenceCost(as.numeric(input$testPeriods), inputParams) 
   # })   
   # 

}

# Run the application 
shinyApp(ui = ui, server = server)

