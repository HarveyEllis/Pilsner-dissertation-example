// (C) P.J. Dodd (2018): Distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include "eventlogic.h"

const double largeNumber = 1000000000; //there will be no one who lives longer 100, but 200 just to be safe

// Other custom functions that have been included - these mirror the simul8 VL sections
// Event times //
double drawVTETime(double now, person& P, NumericVector& Z);
double drawHaemorrhageTime(double now, person& P, NumericVector& Z);
double drawNaturalDeathTime(double now, person& P, NumericVector& Z);
// Update characteristics //
void updatePatientCharacteristics(double now, person& P,NumericVector& Z);
// Outcome distributions
int VTEOutcomesDist(double now, person& P, NumericVector& Z);
int HaemOutcomesDist(double now, person& P, NumericVector& Z);
// Discounting //
double OneoffDiscountFactor(double Time,double discountRate);
double ContinuousDiscountFactor(double startTime, double endTime, double discountRate);
// Debugging //
void debugPatient(person& P);

int eventdefn( person& P,       // person to be acted on
               double& now,     // reference to current time
               event & ev,      // event to be enacted
               eventq& EQ,      // event queue to add consequences to
               NumericVector& Z// ,          // input parameters
                 // gsl_rng *r
){
  now = ev.t;                   // advance time
  if(P.D["alive"]==1){          // RIP - we don't want to act on patients who have already died
    event todo = {.t=now, .who=P.who, .what=""}; // event to be determined
    
    // --------- DEFINE LOGIC BELOW HERE (& possibly the "alive" guard): ---------
    // NB users completely responsible for any eligibility/guard logic & consistency between state/input names!!
    if(ev.what=="initialize"){
      todo.what = "VTE_recurrence"; todo.t = drawVTETime(now,P,Z);
      P.C["timeofVTErecurrence"] = todo.t;
      EQ.push(todo);
      
      todo.what = "haemorrhage"; todo.t = drawHaemorrhageTime(now,P,Z);
      P.C["timeofhaemorrhage"] = todo.t;
      EQ.push(todo);
      
      todo.what="death_natural"; todo.t = drawNaturalDeathTime(now,P,Z);  // time-to-death
      P.C["timeofnaturaldeath"] = todo.t;
      EQ.push(todo);
    }
    if(ev.what=="VTE_recurrence"){
      //update patient characteristics
      updatePatientCharacteristics(now,P,Z);
      
      bool trtRestarted;
      //check diagnosis and put on treatment
      if((P.D["previous_haemorrhage"]!=1) && (R::runif(0,1)<=Z["DiagRecurrence"])){
        if(P.C["trt_subs_endtime"]<now){
          P.C["trt_1Q_endtime"]=now+90.0/365.0;
          trtRestarted = true;
          
        }
        P.C["trt_subs_endtime"]=largeNumber;  
      } else {
        P.D["VTE_highrisk"]=1; //this is simply a flag for overwriting purposes when VTE gets redrawn
        trtRestarted = false;
      }
      
      //do outcomes
      int VTEOutcome = VTEOutcomesDist(now,P,Z); //new var as this is used below too
      switch(VTEOutcomesDist(now,P,Z)){
      case 1: //Patient recovers
        P.C["costs"]+= Z["DVTCost"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]); 
        break; 
      case 2: //Patient has PTS
        P.C["costs"]+= Z["PTSCost"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]);
        P.D["previous_PTS"] = 1;
        break; 
      case 3: //Patient has a non-fatal PE
        P.C["costs"]+= Z["PECostNonFatal"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]);
        P.D["previous_PE"] = 1;
        break; 
      case 4: //Patient has a fatal PE
        P.C["costs"]+= Z["PECostFatal"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]); 
        P.D["cause_of_death"]= 2;
        P.C["timeofdeath"] = ev.t;  
        P.D["alive"]=0;
        break; 
      }
      
      //draw other events if the patient survives
      if(VTEOutcome != 4){
        //draw new event and add it
        todo.what = "VTE_recurrence"; todo.t = drawVTETime(now,P,Z);
        P.C["timeofVTErecurrence"] = todo.t;
        EQ.push(todo);
        
        if (trtRestarted){
          todo.what = "haemorrhage"; todo.t = drawHaemorrhageTime(now,P,Z);
          P.C["timeofhaemorrhage"] = todo.t;
          EQ.push(todo);
        }
        
      }
    }
    
    if(ev.what=="haemorrhage"){
      //update
      updatePatientCharacteristics(now,P,Z);
      
      int HaemOutcome = HaemOutcomesDist(now,P,Z); //new var as this is used below too
      switch(HaemOutcome){
      case 1: //fatal
        P.C["costs"]+= Z["HaemCostFatal"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]); 
        P.D["cause_of_death"]=3;
        P.C["timeofdeath"] = ev.t;  
        P.D["alive"]=0;
        break; 
      case 2: //intracranial
        P.C["costs"]+= Z["HaemCostIntracranialImpact"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]);
        P.D["previous_haemorrhage"] = 1;
        break; 
      case 3: //non-intracranial
        P.C["costs"]+= Z["HaemCostNonIntracranial"] * OneoffDiscountFactor(now, Z["DiscCostsOneOff"]);
        P.D["previous_haemorrhage"] = 2;
        break; 
      }
      //cease treatment
      P.C["trt_subs_endtime"]=now;
      P.C["trt_1Q_endtime"] = min(P.C["trt_1Q_endtime"],now);
      
      //draw other events if the patient survives
      if(HaemOutcome != 1){
        if(P.D["VTE_highrisk"]==0){
          todo.what = "VTE_recurrence"; todo.t = drawVTETime(now,P,Z);
          P.C["timeofVTErecurrence"] = todo.t;
          EQ.push(todo);
        }
        
        todo.what = "haemorrhage"; todo.t = drawHaemorrhageTime(now,P,Z);
        P.C["timeofhaemorrhage"] = todo.t;
        EQ.push(todo);
      }
    }
    
    if(ev.what=="death_natural"){
      updatePatientCharacteristics(now,P,Z);
      P.D["cause_of_death"]=1;
      P.C["timeofdeath"] = ev.t;            // recorded as dead
      P.D["alive"] = 0;              // recorded as dead
      
    }
    if(ev.what=="finalize"){
      // no other changes needed: just for updating ages if alive
    }
    // --------- END LOGIC DEFINITION ---------
    
    updateages(P.A,ev.t-P.tle);    // update all age-like things
    P.tle = ev.t;       // doesn't matter if no change; last event time used for aging
    
  }
  return 0;
  
}

// ------------------ Event draws ------------------ //
double drawVTETime(double now, person& P, NumericVector& Z){
  double eventTime;
  double EndTrtTime = P.C["trt_subs_endtime"];
  
  //check for on treatment
  if(now<=EndTrtTime){ //if patient is on trt then draw trt time
    eventTime = now + R::rexp(1/Z["VTERiskTrtLambda"]);
    
    if(eventTime>EndTrtTime){
      eventTime = EndTrtTime+R::rexp(1/Z["VTERiskUnTrtLambda"]);
    }
  } else { //patient not on treatment
    eventTime = now + R::rexp(1/Z["VTERiskFirst6moLambda"]);
    if(eventTime>now + 180.0/365.0){ //if it is more than 6 months, then redraw and set to normal risk
      eventTime = now + 180.0/365.0+R::rexp(1/Z["VTERiskUnTrtLambda"]);
    }
  }
  return now + eventTime;
};

double drawHaemorrhageTime(double now, person& P, NumericVector& Z){
  //check first three months
  double eventTime = now + R::rexp(1/Z["HaemRiskInit3moLambda"]);
  
  if(eventTime > P.C["trt_1Q_endtime"]){
    //Draw the first event
    double drawnTime = R::rexp(1/Z["HaemRiskLessThan40Lambda"]);
    int i = 30; //starting age
    
    double accumulatedTime = 0;
    double maxNoYearsInBand = 0;
    while(i<=70){
      
      
      if(i==70){
        maxNoYearsInBand = largeNumber;
      } else{
        maxNoYearsInBand = max((double)i+10.0 - max(P.C["age"],(double)i),0.0); //rem to put 0.0 to force type to double
      }
      
      //'MIN of the random draw and maximum number of years is the number of years added to the accumulated years
      accumulatedTime += min(drawnTime,maxNoYearsInBand); //if there can be 0 years in the band then 0 years are added
      
      //we then check if we need a redraw because the drawn time is beyond the max no of years in the band (it always will be if the band is 0)
      if(drawnTime > maxNoYearsInBand){
        switch(i){
        case 30: drawnTime = R::rexp(1/Z["HaemRisk40to49Lambda"]); break;
        case 40: drawnTime = R::rexp(1/Z["HaemRisk50to59Lambda"]); break;
        case 50: drawnTime = R::rexp(1/Z["HaemRisk60to69Lambda"]); break;
        case 60: drawnTime = R::rexp(1/Z["HaemRiskMoreThan70Lambda"]); break;
        }
      } else{
        break;
      }
      i += 10;
    }
    //set the event time
    eventTime = now + accumulatedTime;
  }
  
  // //There is no risk from haemorrhage after the patient comes off treatment
  if(eventTime>P.C["trt_subs_endtime"]){eventTime = largeNumber;}
  return eventTime;
};

double drawNaturalDeathTime(double now, person& P, NumericVector& Z){
  string str;
  int deathAge = 0;
  //should always work since we only draw this time at the start
  
  int i = (int) P.C["age"];
  while(i<=100){
    str = "qx"+to_string(i);
    if(R::runif(0,1)<= Z[str]){
      deathAge = i;
      break;
    } else if(i == 100){
      deathAge = 100;
      break;
    }
    i++;
  }
  return (double) now + deathAge-P.C["age"]+R::runif(0,1);
}

// ------------------ Updating patient characteristics ------------------ //
// Sub routine to update patient characteristics
void updatePatientCharacteristics(double now, person& P,NumericVector& Z){
  //apply initial period costs
  double starttime = max(min(now, P.C["trt_starttime"]),P.tle);
  double endtime = max(min(now, P.C["trt_1Q_endtime"]),P.tle);
  
  P.C["costs"] += Z["TrtCostImplementation"]*(365.0/90.0)* ContinuousDiscountFactor(starttime,endtime,Z["DiscCostsContinuous"]);
  
  //apply subsequent treatment costs
  starttime = max(min(now, P.C["trt_1Q_endtime"]),P.tle);
  endtime = max(min(now, P.C["trt_subs_endtime"]),P.tle);
  
  P.C["costs"] += Z["TrtCostMaintenance"]* ContinuousDiscountFactor(starttime,endtime,Z["DiscCostsContinuous"]);
  
  //update other costs
  if(P.D["previous_haemorrhage"]==1){
    P.C["costs"] += Z["HaemCostIntracranialMaintenance"]* ContinuousDiscountFactor(P.tle,now,Z["DiscCostsContinuous"]);
  }
  
  //turn off high risk period if applicable
  if(now - P.tle > 180.0/365.0){
    P.D["VTE_highrisk"]=0;
  }
  
  //update utilities since last event
  int i = 0;
  while((i<=70) && (i<=now)){ //we do this based around time instead of age!
    double bottomAgeBand = max(min(now,(double)i),P.tle);
    double topAgeBand = max(min(now,(double)i+10),P.tle);
    double endTrtMidPoint = max(min(P.C["trt_subs_endtime"],topAgeBand),bottomAgeBand);
    double utility = 0;
    
    //get utility based on the model time and the age the patient entered the model
    int switcher = min(i+P.D["age_group"],70); //so that we can't go over
    switch(switcher){
    case 30:
      utility = Z["UtilBase30to39"];
      break;
    case 40:
      utility = Z["UtilBase40to49"];
      break;  
    case 50:
      utility = Z["UtilBase50to59"];
      break;
    case 60:
      utility = Z["UtilBase60to69"];
      break;
    case 70:
      utility = Z["UtilBase70plus"];
      break;
    }
    
    if(P.D["previous_PE"]==1){ //previous PE
      utility *= Z["UtilPENonFatal"];
    }
    if(P.D["previous_haemorrhage"]==1){ //previous intracranial
      utility *= Z["UtilHaemIntracranial"];
    } else if (P.D["previous_haemorrhage"]==2){ //previous nonintracranial
      utility *= Z["UtilHaemNonIntracranial"];
    }
    if(P.D["previous_PTS"]==1){
      utility *= Z["UtilPTS"];
    }
    
    //apply non treated utility
    utility *= Z["UtilOffTrt"];
    //P.C[""]
    P.C["QALYs"] += utility * ContinuousDiscountFactor(endTrtMidPoint,topAgeBand,Z["DiscQALYsContinuous"]);
    //apply treated utility
    utility *= Z["UtilOnTrt"];
    P.C["QALYs"] += utility * ContinuousDiscountFactor(bottomAgeBand,endTrtMidPoint,Z["DiscQALYsContinuous"]);
    i += 10;
  }
  //update ages
  P.C["age"]= (double)P.D["age_group"] + now;
};

// ------------------ Event Outcomes ------------------ //
// Draws for the outcomes of the events
int VTEOutcomesDist(double now, person& P, NumericVector& Z){
  IntegerVector outcomeframe = IntegerVector::create(1, 2, 3, 4);
  NumericVector probs;
  
  if(P.D["previous_PE"]==1){
    if (now<=P.C["trt_subs_endtime"]){
      probs = {Z["VTEOutcomeTrtPEtoDVT"],Z["VTEOutcomeTrtPEtoPTS"],Z["VTEOutcomeTrtPEtoPE"],Z["VTEOutcomeTrtPEtoDead"]};
    } else{
      probs = {Z["VTEOutcomeUnTrtPEtoDVT"],Z["VTEOutcomeUnTrtPEtoPTS"],Z["VTEOutcomeUnTrtPEtoPE"],Z["VTEOutcomeUnTrtPEtoDead"]};
    }
  } else {
    if (now<=P.C["trt_subs_endtime"]){
      probs = {Z["VTEOutcomeTrtDVTtoDVT"],Z["VTEOutcomeTrtDVTtoPTS"],Z["VTEOutcomeTrtDVTtoPE"],Z["VTEOutcomeTrtDVTtoDead"]};
    } else{
      probs = {Z["VTEOutcomeUnTrtDVTtoDVT"],Z["VTEOutcomeUnTrtDVTtoPTS"],Z["VTEOutcomeUnTrtDVTtoPE"],Z["VTEOutcomeUnTrtDVTtoDead"]};
    }
  }
  return sample(outcomeframe,1,true,probs)[0]; //take the first element of the integer vector to return as an int.
}

int HaemOutcomesDist(double now, person& P, NumericVector& Z){
  IntegerVector outcomeframe = IntegerVector::create(1, 2, 3);
  NumericVector probs;
  if(now<=P.C["trt_1Q_endtime"]){
    probs = {Z["HaemOutcomeInit3moDead"],Z["HaemOutcomeInit3moIntracranial"],Z["HaemOutcomeInit3moNonintracranial"]};
  } else {
    probs = {Z["HaemOutcomeSubsqDead"],Z["HaemOutcomeSubsqIntracranial"],Z["HaemOutcomeSubsqNonintracranial"]};
  }
  return sample(outcomeframe,1,true,probs)[0]; //take the first element of the integer vector to return as an int.
};

// ------------------ Discounting ------------------ //
double ContinuousDiscountFactor(double startTime, double endTime, double discountRate){
  return (exp(endTime*-discountRate)-exp(startTime*-discountRate))/-discountRate;
}

double OneoffDiscountFactor(double Time,double discountRate){
  return 1/(pow(1+discountRate,Time));
}

// ------------------ Debugging ------------------ //
void debugPatient(person& P){
  auto debugList = NumericVector::create(
    _["Time in model"] = P.A["timeinmodel"],
    _["age"] = P.C["age"],
    _["Treatment start time"] = P.C["trt_starttime"],
    _["Treatment end time 1Q"] = P.C["trt_1Q_endtime"],
    _["Treatment start time subsequently"] = P.C["trt_subs_endtime"],
    _["Costs"] = P.C["costs"],              
    _["QALYs"] = P.C["QALYs"],    
    _["Time of death"] = P.C["timeofdeath"],
    _["Time of VTE recurrence"] = P.C["timeofVTErecurrence"],
    _["Time of haemorrhage"] = P.C["timeofhaemorrhage"],
    _["Cause of death"] = P.D["cause_of_death"],
    _["Sex"] = P.D["sex"],
    _["Previous PE"] = P.D["previous_PE"],
    _["Previous haemorrhage"] = P.D["previous_haemorrhage"],
    _["Previous PTS"] = P.D["previous_PTS"],
    _["Thrombophilia type"] = P.D["thrombophilia_type"],
    _["High risk VTE period"] = P.D["VTE_highrisk"],
    _["Alive"] = P.D["alive"]);
  Rf_PrintValue(debugList);
}