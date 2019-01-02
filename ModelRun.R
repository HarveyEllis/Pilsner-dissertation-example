#### Load in parameters, populations and compile
## Library and compile
library(Rcpp) # assuming in dir above src
library(tidyverse)
library(magrittr)
library("PILSNERSim")
library(foreach)
library(doParallel)

#### Compiles ####
#Note you must set up a local package first, using a command such as RcppArmadillo.package.skeleton( "PILSNERSim" ) 
#Then install using the following command
#install.packages("PILSNERSim", repos=NULL, type="source")
#sourceCpp('src/pilsner.cpp') #compiles: only needed once per session

#### Settings and parameters #### 
## Model modes
ParameterMode <- "updated" # "original" or "updated"
RunMode<- "probabilistic" # "deterministic" or "probabilistic"

if(ParameterMode == "original") {
  Populations<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Model populations original.csv",check.names = FALSE)
  Parameters<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Parameters Original.csv", 
                       check.names = FALSE)
  Lifetables<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Lifetables 2003-2005.csv",check.names = TRUE, fileEncoding="UTF-8-BOM")
} else if (ParameterMode == "updated") {
  Populations<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Model populations updated.csv",check.names = FALSE)
  Parameters<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Parameters Updated.csv", 
                       check.names = FALSE)
  Lifetables<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Lifetables 2015-2017.csv",check.names = TRUE, fileEncoding="UTF-8-BOM")
}

if(RunMode == "probabilistic") {
  Populations<-read.csv("C:/Users/Harvey/Files/University/Thesis/Model params and populations/Model populations PSA.csv",check.names = FALSE)
}

MaleLifeTable<-Lifetables[["MaleQx"]]%>%set_names(paste0("qx",Lifetables[["Age"]]))
FemaleLifeTable<-Lifetables[["FemaleQx"]]%>%set_names(paste0("qx",Lifetables[["Age"]]))

ModelLargeNumber<-1000000000
NumPopRuns<-nrow(Populations)
NumPSARuns<-2500
NoPatients <-2500 #cohort size

DiscountCostsPerYear<-0.035
DiscountQALYsPerYear<-0.035
ContinuousDiscountRateCosts<- -log(1/(1+DiscountCostsPerYear))
ContinuousDiscountRateQALYs<- -log(1/(1+DiscountQALYsPerYear))

#### Setup parameters into the correct form  ####
# parameters for risk of recurrence of VTE
LoadPopulationAndParameters<-function(PopID,RunID,RunMode){
  if(RunMode=="probabilistic"){RunID<-RunID+1} ## create a quick offset for the parameters - deterministic is in the first row, and PSA is after that
  
  # this is not really used as the other matrices suffice
  AM <<- matrix(0.00000000000000000000000001,ncol=1,nrow=NoPatients, dimnames = list(c(),c("timeinmodel"))); 
  
  CM <<- matrix(c(Populations[PopID,"Age"],
                  0,
                  90/365,
                  switch (EXPR = Populations[PopID,"Duration of treatment"],
                          '1' = ifelse(Populations[PopID,"Index Event"]==1, 90/365,180/365),
                          '2' = 10,
                          '3' = 20,
                          '4' = ModelLargeNumber
                  ),
                  0,
                  0, 
                  0,   
                  0,
                  0,
                  0),
                ncol=10,nrow=NoPatients,byrow = TRUE, 
                dimnames = list(c(),c("age",  ##this is incremented
                                      "trt_starttime",
                                      "trt_1Q_endtime",
                                      "trt_subs_endtime",
                                      "costs",
                                      "QALYs",
                                      "timeofdeath",
                                      "timeofnaturaldeath",
                                      "timeofVTErecurrence",
                                      "timeofhaemorrhage")))
  DM <<- matrix(c(
    Populations[PopID,"Age"],
    0,                                            
    Populations[PopID,"Sex"],
    ifelse(Populations[PopID,"Index Event"]==1,0,1),
    0,
    0,
    Populations[PopID,"Thrombotype"],
    0,
    1),
    ncol=9,nrow=NoPatients,byrow = TRUE,
    dimnames= list(c(),c("age_group", ##this is just for holding what age patients started in
                         "cause_of_death",
                         "sex",
                         "previous_PE",
                         "previous_haemorrhage",
                         "previous_PTS",
                         "thrombophilia_type",
                         "VTE_highrisk",
                         "alive")))
  
  ## Calculations for VTE risk of recurrence
  RelativeRiskVTEThrombotype<-switch (EXPR = Populations[PopID,"Thrombotype"], 
                                      '1' = 1,
                                      '2' = Parameters[RunID, "Relative risk of recurrence Lupus anticoagulant"],
                                      '3' = Parameters[RunID,"Relative risk of recurrence FVL heterozygous and PTG20210A heterozygous"],
                                      '4' = Parameters[RunID,"Relative risk of recurrence Anticardiolipin antibody"],
                                      '5' = Parameters[RunID,"Relative risk of recurrence FVL homozygous"],
                                      '6' = Parameters[RunID,"Relative risk of recurrence PTG20210A heterozygous"],
                                      '7' = Parameters[RunID,"Relative risk of recurrence AT deficiency PC deficiency or PS deficiency or lupus-like anticoagulants"],
                                      '8' = Parameters[RunID,"Relative risk of recurrence FVL heterozygous"])
  
  VTERecurrenceBaselineLambda<-ifelse(Populations[PopID, "Sex"]==1,
                                      Parameters[RunID,"Risk of recurrent VTE in men lambda"],
                                      Parameters[RunID,"Risk of recurrent VTE in women lambda"])
  VTERecurrenceLambdaUntreated<- VTERecurrenceBaselineLambda*RelativeRiskVTEThrombotype # simple rate ratio
  VTERecurrenceLambdaTreated<-VTERecurrenceLambdaUntreated*(1-Parameters[RunID,"Efficacy of warfarin on VTE recurrence"])
  
  HaemorrhageHazardRatio<-ifelse(Populations[PopID,"Treatment"]==1, 
                                 1, 
                                 Parameters[RunID,"Efficacy of NOACs vs warfarin on bleeding"])
  
  ## Other variables that change based on the population
  TrtCostImplementation<-ifelse(Populations[PopID,"Treatment"]==1, 
                                Parameters[RunID,"Implementation of warfarin therapy"], 
                                Parameters[RunID,"Implementation cost of NOAC therapy"])
  
  TrtCostMaintenance<-ifelse(Populations[PopID,"Treatment"]==1, 
                             Parameters[RunID,"Maintenance of warfarin therapy yearly"], 
                             Parameters[RunID,"Maintenance of NOAC therapy yearly"])
  
  TrtUtil<-ifelse(Populations[PopID,"Treatment"]==1, 
                  Parameters[RunID,"Post VTE Receiving warfarin"], 
                  Parameters[RunID,"Post VTE Receiving NOACs"])
  
  ## Life tables
  CurrentLifetable<-if(Populations[PopID, "Sex"]==1){MaleLifeTable} else {FemaleLifeTable}
  
  ## Add all other variables that don't depend on the population to the dataframe - nothing fancy, just type them out for clarity
  PP<<-c("VTERiskTrtLambda" = VTERecurrenceLambdaTreated,
         "VTERiskUnTrtLambda" = VTERecurrenceLambdaUntreated,
         "VTERiskFirst6moLambda" = Parameters[RunID,"Risk of recurrence in first 6 mo lambda"],
         "VTEOutcomeTrtDVTtoDVT" = Parameters[RunID,"Treated DVTs that result in resolving DVT at next VTE"],
         "VTEOutcomeTrtDVTtoPTS" = Parameters[RunID,"Treated DVTs that result in PTS at next VTE"],
         "VTEOutcomeTrtDVTtoPE" = Parameters[RunID, "Treated DVTs that result in non-fatal PE at next VTE"],
         "VTEOutcomeTrtDVTtoDead" = Parameters[RunID, "Treated DVTs that result in fatal PE at next VTE"],
         "VTEOutcomeUnTrtDVTtoDVT" = Parameters[RunID,"Untreated DVTs that result in resolving DVT at next VTE"],
         "VTEOutcomeUnTrtDVTtoPTS" = Parameters[RunID,"Untreated DVTs that result in PTS at next VTE"],
         "VTEOutcomeUnTrtDVTtoPE" = Parameters[RunID, "Untreated DVTs that result in non-fatal PE at next VTE"],
         "VTEOutcomeUnTrtDVTtoDead" = Parameters[RunID, "Untreated DVTs that result in fatal PE at next VTE"],
         "VTEOutcomeTrtPEtoDVT" = Parameters[RunID,"Treated PEs that result in resolving DVT at next VTE"],
         "VTEOutcomeTrtPEtoPTS" = Parameters[RunID,"Treated PEs that result in PTS at next VTE"],
         "VTEOutcomeTrtPEtoPE" = Parameters[RunID, "Treated PEs that result in non-fatal PE at next VTE"],
         "VTEOutcomeTrtPEtoDead" = Parameters[RunID, "Treated PEs that result in fatal PE at next VTE"],
         "VTEOutcomeUnTrtPEtoDVT" = Parameters[RunID,"Untreated PEs that result in resolving DVT at next VTE"],
         "VTEOutcomeUnTrtPEtoPTS" = Parameters[RunID,"Untreated PEs that result in PTS at next VTE"],
         "VTEOutcomeUnTrtPEtoPE" = Parameters[RunID, "Untreated PEs that result in non-fatal PE at next VTE"],
         "VTEOutcomeUnTrtPEtoDead" = Parameters[RunID, "Untreated PEs that result in fatal PE at next VTE"],
         "HaemRiskInit3moLambda" = Parameters[RunID, "Probability of haemorrhage in initial 3 months of treatment lambda"] * HaemorrhageHazardRatio,
         "HaemRiskLessThan40Lambda" = Parameters[RunID, "Probability of haemorrhage between the ages of 30 and 39 lambda"] * HaemorrhageHazardRatio,
         "HaemRisk40to49Lambda" = Parameters[RunID, "Probability of haemorrhage between the ages of 40 and 49 lambda"] * HaemorrhageHazardRatio,
         "HaemRisk50to59Lambda" = Parameters[RunID, "Probability of haemorrhage between the ages of 50 and 59 lambda"] * HaemorrhageHazardRatio,
         "HaemRisk60to69Lambda" = Parameters[RunID, "Probability of haemorrhage between the ages of 60 and 69 lambda"] * HaemorrhageHazardRatio,
         "HaemRiskMoreThan70Lambda" =Parameters[RunID, "Probability of haemorrhage at age 70 and above"] * HaemorrhageHazardRatio,
         "HaemOutcomeInit3moDead" = Parameters[RunID, "Haemorrhage in initial 3 months that are fatal"],
         "HaemOutcomeInit3moIntracranial" = Parameters[RunID, "Haemorrhage in initial 3 months that are non-fatal intracranial"],
         "HaemOutcomeInit3moNonintracranial" = Parameters[RunID, "Haemorrhage in initial 3 months that are non-fatal non-intracranial"],
         "HaemOutcomeSubsqDead" = Parameters[RunID, "Haemorrhage subsequently that are fatal"],
         "HaemOutcomeSubsqIntracranial" = Parameters[RunID, "Haemorrhage subsequently that are non-fatal intracranial"],
         "HaemOutcomeSubsqNonintracranial" = Parameters[RunID, "Haemorrhage subsequently that are non-fatal non-intracranial"],
         "DiagRecurrence" = Parameters[RunID, "Diagnosis of DVT"],
         "TrtCostImplementation" = TrtCostImplementation,
         "TrtCostMaintenance" = TrtCostMaintenance,
         "DVTCost" = Parameters[RunID, "DVT Cost of further resolving"],
         "PTSCost" = Parameters[RunID, "PTS Cost of treating"],
         "PECostFatal" = Parameters[RunID, "Treating a fatal PE overall cost"],
         "PECostNonFatal" = Parameters[RunID, "Treating a non-fatal PE overall cost"],
         "HaemCostFatal" = Parameters[RunID, "Treating a fatal haemorrhage Standard unit cost"],
         "HaemCostIntracranialImpact" = Parameters[RunID, "Treating a non-fatal intracranial haemorrhage Impact cost"],
         "HaemCostIntracranialMaintenance" = Parameters[RunID, "Treating a non-fatal intracranial haemorrhage Maintenance cost yearly"],
         "HaemCostNonIntracranial" = Parameters[RunID, "Treating a non-fatal non-intracranial haemorrage overall cost"],
         "UtilOnTrt" =  TrtUtil,
         "UtilOffTrt" = Parameters[RunID, "Post VTE Not Receiving warfarin"],
         "UtilPTS" = Parameters[RunID, "Post-thrombotic syndrome"],
         "UtilPENonFatal" = Parameters[RunID, "Non-fatal pulmonary embolism"],
         "UtilHaemIntracranial" = Parameters[RunID, "Non-fatal intracranial haemorrhage"],
         "UtilHaemNonIntracranial" = Parameters[RunID, "Non-fatal non-intracranial haemorrhage"],
         "UtilBase30to39" = Parameters[RunID, "Baseline 30 to 39"],
         "UtilBase40to49" = Parameters[RunID, "Baseline 40 to 49"],
         "UtilBase50to59" = Parameters[RunID, "Baseline 50 to 59"],
         "UtilBase60to69" = Parameters[RunID, "Baseline 60 to 69"],
         "UtilBase70plus" = Parameters[RunID, "Baseline 70 plus"],
         CurrentLifetable,
         "DiscCostsContinuous" = ContinuousDiscountRateCosts,
         "DiscCostsOneOff" = DiscountCostsPerYear,
         "DiscQALYsContinuous" = ContinuousDiscountRateQALYs,
         "DiscQALYsOneOff" = DiscountQALYsPerYear,
         "ModelLargeNumber" = ModelLargeNumber)
  Params<<-rep(list(PP),each=NoPatients)
}

#### Main Loop - For each implementation ####
registerDoParallel(makeCluster(8))  
opts <- list(chunkSize=10)
#N.B: Don't try and use side-effects from a foreach parallel loop - you have to return everything!!
start_time<-proc.time()[[3]]
Results<-foreach(RunID=1:ifelse(RunMode=="deterministic",1,NumPSARuns), .combine='rbind', .options.nws=opts) %:% 
  foreach (PopID=1:NumPopRuns,.combine = rbind, .packages = c("PILSNERSim","dplyr","magrittr")) %dopar% {
    LoadPopulationAndParameters(PopID, RunID, RunMode)    
    set.seed(123)
    PILSNERSim::simulator(AM,CM,DM,ModelLargeNumber,Params,recording=FALSE)    
    data.frame("RunID" = RunID, "PopID" =  PopID, "costs"=sum(CM[,"costs"]), "QALYs" = sum(CM[,"QALYs"]))
  }
end_time<-proc.time()[[3]]-start_time
stopImplicitCluster() #remember to stop the cluster to avoid leaks!!

#### Run time testing ####
registerDoParallel(makeCluster(8))  
opts <- list(chunkSize=10)
perfresults<-data.frame()
for(i in 1:5){
  start_time<-proc.time()[[3]]
  Results<-foreach(RunID=1:1, .combine = rbind, .options.nws=opts) %:% 
    foreach (PopID=1:640,.combine = rbind, .packages = c("PILSNERSim","dplyr","magrittr")) %dopar% {
      LoadPopulationAndParameters(PopID, RunID, "deterministic")    
      set.seed(123)
      PILSNERSim::simulator(AM,CM,DM,ModelLargeNumber,Params,recording=FALSE)    
      data.frame("RunID" = RunID, "PopID" =  PopID, "costs"=sum(CM[,"costs"]), "QALYs" = sum(CM[,"QALYs"]))
    }
  end_time<-proc.time()[[3]]-start_time
  perfresults%<>%bind_rows(data.frame(end_time))
}
stopImplicitCluster() #remember to stop the cluster to avoid leaks!!

#### Print out results ####
write.csv(x = Results, file = "C:/Users/Harvey/Files/University/Thesis/Model results/R Output results.CSV")
write.csv(x = perfresults, file = "C:/Users/Harvey/Files/University/Thesis/Model results/R Output perf results.CSV")