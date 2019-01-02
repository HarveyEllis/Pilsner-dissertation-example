// (C) P.J. Dodd (2018): Distributed under CC BY 4.0 license https://creativecommons.org/licenses/by/4.0/
#ifndef EVENTDEFN               /* guard */
#define EVENTDEFN

// [[Rcpp::plugins("cpp11")]]
#include <iostream>
#include <stdlib.h>
#include <queue>
#include <fstream>
#include <ctime>
#include <math.h>
#include <vector>
#include <map>
//#include <Rcpp.h> - changed by HE 
#include "RcppArmadillo.h"


/* /\* below to get GSL working *\/ */
/* #include <RcppGSL.h> */
/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_randist.h> */

using namespace std;
using namespace Rcpp;

/* person class */
class person{
 public:
  unsigned int who;             // population id
  double tle;                   // time of last event
  map< string, double> A;       // ages: age, other times-since-x (that advance with time)
  map< string, double> C;       // continuous states
  map< string, int> D;          // integer states
  // needs constructor
  person(unsigned int id,
         vector< string > AL, vector< double > AV,
         vector< string > CL, vector< double > CV,
         vector< string > DL, vector< int > DV);
};


// event structure
struct event{
  double t;
  unsigned int who;
  string what;
};

// comparator to make the priority queue work
struct LTbyTime{
  bool operator()(const event& lhs, const event& rhs) const {
    return lhs.t > rhs.t;       //care with direction!
  }
};

// for event queue
typedef priority_queue< event, vector<event>, LTbyTime> eventq;

// generic approach for converting maps to vectors
template <typename M, typename V>
  void MapToVec( const  M & m, V & v ) {
  for( typename M::const_iterator it = m.begin(); it != m.end(); ++it ) {
    v.push_back( it->second );
  }
}


// results class here
class results{
public:
  bool recording;

  // --- data ------
  vector < double > tz;
  vector< unsigned int> whoz;
  vector< string > whatz;
  vector< vector< double> > Az;
  vector< vector< int > > Dz;
  vector< vector< double> > Cz;

  // ---- functions ----

  // function to capture event from event and person
  int record(event& todo, person& P){
    // need to convert the map to vector
    tz.push_back(todo.t);
    whoz.push_back(todo.who);
    whatz.push_back(todo.what);
    vector< double > dvec;
    vector< int > ivec;
    MapToVec< map< string, double>, vector< double > > (P.A,dvec); // convert
    Az.push_back(dvec);
    MapToVec< map< string, int>, vector< int > > (P.D,ivec); // convert
    Dz.push_back(ivec);
    dvec.clear();                                                  // empty out
    MapToVec< map< string, double>, vector< double > > (P.C,dvec); // convert
    Cz.push_back(dvec);
    return 0;
  }

  DataFrame getresult(){        //
    // depends on res.recording
    DataFrame ans;
    if(recording){
      /* need to onvert tz to Numeric vector */
      NumericVector tzn; tzn.assign(tz.begin(),tz.end());
      StringVector wtzn; wtzn.assign(whatz.begin(),whatz.end());
      IntegerVector wozn; wozn.assign(whoz.begin(),whoz.end());
      ans = DataFrame::create(Named("time")=tzn,Named("event")=wtzn,Named("who")=wozn);
      /* TODO think about adding on state data for debugging, reshaping bit of pain */
    } else
      ans = DataFrame::create(Named("recording")=0,Named("turned")=0,Named("off")=0);
    return ans;
  }

}; /* will need function defns in pilsner */

// =========== function prototypes ============
int eventdoer( eventq& EQ, vector< person >& cohort, results& res,
               double& now, List& parmzl );

int updateages( map< string, double >& mymap, double dt );

/* actual event defn declaration */
int eventdefn( person& P,
               double& now,
               event & ev,
               eventq& EQ,
               NumericVector& Z);/* , */
               /* gsl_rng *r); */
#endif
