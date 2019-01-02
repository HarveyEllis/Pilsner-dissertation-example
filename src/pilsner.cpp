// (C) P.J. Dodd (2018): Distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
#include "eventlogic.h"         // includes other class definitions etc
// [[Rcpp::plugins("cpp11")]]

// =========== class definitions ===============

// person class declared in eventlogic

// constructor
person::person(unsigned int id,
               vector< string >  AL, vector< double > AV,
               vector< string >  CL, vector< double > CV,
               vector< string >  DL, vector< int > DV){
  who=id;
  for(unsigned int i = 0; i < AL.size(); ++i ){
    A[AL.at(i)] = AV.at(i);
  }
  for(unsigned int i = 0; i < DL.size(); ++i ){
    D[DL.at(i)] = DV.at(i);
  }
  for(unsigned int i = 0; i < CL.size(); ++i ){
    C[CL.at(i)] = CV.at(i);
  }
}

//=========== MAIN SIMULATOR ======================
//[[depends(RcppGSL)]] to be altered if attempting GSL PRNGs
// [[Rcpp::export]]
DataFrame simulator( NumericMatrix& A, NumericMatrix& C,
               NumericMatrix& D,// ages, continuous states, discrete
               double endtime,                      // when to stop
               List& parmzl,                 // parameters for model
               bool recording = false // whether to record all events
               ){

  // --------- declarations etc ---------
  if( (A.nrow() != D.nrow()) | (A.nrow() != C.nrow()) | (C.nrow() != D.nrow())) stop("A, D, C must have the same number of rows!");
  unsigned int N = (unsigned int) C.nrow();               // cohort size
  double stime(0);                   // start time
  eventq EQ;                         // event queue
  results res;                       // results store
  //clock_t start;                     // timer             //commented out by HE
  res.recording = recording;
  vector< person > cohort;           // all the people
  unsigned int ecount(0);            // event count
  CharacterVector Anmz0, Dnmz0, Cnmz0;
  vector< string > Anmz, Dnmz, Cnmz;
  //start = clock();              // time start             //commented out by HE
  Anmz0 = colnames(A); Dnmz0 = colnames(D); Cnmz0 = colnames(C);
  for(int i = 0; i < Anmz0.size(); ++i) Anmz.push_back(string(Anmz0[i]));
  for(int i = 0; i < Cnmz0.size(); ++i) Cnmz.push_back(string(Cnmz0[i]));
  for(int i = 0; i < Dnmz0.size(); ++i) Dnmz.push_back(string(Dnmz0[i]));

  // gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

  // --------- create cohort ---------
  vector< double > AV; vector< int > DV; vector< double > CV;
  for( unsigned int i = 0; i<N; ++i ){ // rewrite with NumericMatrix::Row    row = m( 1 , _ ); ?
    int ii = (int) i;
    AV.clear(); DV.clear(); CV.clear();
    for( int j = 0; j < A.ncol(); ++j ) AV.push_back(A(ii,j));
    for( int j = 0; j < D.ncol(); ++j ) DV.push_back(round(D(ii,j))); // TODO integer issue
    for( int j = 0; j < C.ncol(); ++j ) CV.push_back(C(ii,j));
    cohort.push_back( person(i,Anmz,AV,Cnmz,CV,Dnmz,DV) ); // input data to make
  }

  // --------- initialize ---------
  event initialize = {.t=0.0, .who=0, .what="initialize"};
  event finalize = {.t=endtime, .who=0, .what="finalize"};
  for( unsigned int i = 0; i<cohort.size(); ++i ){
    initialize.who = i; finalize.who = i;
    // make sure everyone has an initializ/finalization step queued (may/may not be defined)
    EQ.push(initialize);
    EQ.push(finalize);
  }

  // --------- simulation loop ---------
  while( stime <= endtime && EQ.size()>0 ){
    ecount++;
    eventdoer(EQ,cohort,res,stime,parmzl);
  }

  // --------- output ---------
  // extract final state into matrices
  for( unsigned int i = 0; i<N; ++i ){
    int ii = (int) i;
    for( int j = 0; j < A.ncol(); ++j ) A(ii,j) = cohort.at(i).A[ Anmz[ (unsigned int)j ] ];
    for( int j = 0; j < C.ncol(); ++j ) C(ii,j) = cohort.at(i).C[ Cnmz[ (unsigned int)j ] ];
    for( int j = 0; j < D.ncol(); ++j ) D(ii,j) = (double)cohort.at(i).D[ Dnmz[ (unsigned int)j ] ];
  }

  // gsl_rng_free(r);

  //Commented out by HE
  //Rcout << ecount << " events took place...in ";
  //Rcout << (clock() - start)  / (double) CLOCKS_PER_SEC;
  //Rcout << " seconds."<< endl; //
  return res.getresult();
}


//=========== function definitions ======================

// for actually doing events
int eventdoer( eventq& EQ,
               vector< person >& cohort,
               results& res,
               double& now,
               List& parmzl// ,
               // gsl_rng *r
               ){
  event todo = EQ.top();        // get work
  // cout << todo.what << " --- " << todo.t << endl;
  EQ.pop();                     // remove from list
  NumericVector parmz;
  int whom = (int) todo.who;
  if(parmzl.size()==1) whom = 0;
  parmz = parmzl[whom];
  eventdefn( cohort.at(todo.who), now, todo, EQ, parmz ); // carry out work
  if(res.recording) res.record(todo,cohort.at(todo.who)); // record event
  return 0;
}

// advance all ages function for maps
int updateages( map< string, double >& mymap, double dt ){
  map<string,double>::iterator it;
  for( it = mymap.begin(); it != mymap.end(); ++it){
    if( it->second >= 0){         // will use negative to initialize ages
      it->second += dt; //value
    }
  }
  return 0;
}


