#include "observables.h"
#include <cmath>

observables::observables(read_music_output_files* armof) : rmof(armof){
  event_arena = rmof->get_event_arena();
  rand = new random_gen();
  
  // print to check whether reading is perfect or not ???
  /*
  for(int ii=0; ii<rmof->get_music_pit_bins(); ii++){
    double pt = rmof->get_pt_val_of_bin(ii) ; 
    double v4cos = rmof->get_event(1)->get_pt_differential_vn(4, 0, ii);
    double v4sin = rmof->get_event(1)->get_pt_differential_vn(4, 1, ii);
    std::cout << ii << "  " << pt << "  " << v4cos << "  " << v4sin << std::endl ; 
  }
  */
}


std::vector<int> observables::get_an_event_ensemble(){
  std::vector<int> event_ID_ens;
  for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    int evID =  rmof->get_total_music_events() * rand->rand_uniform() ;
    event_ID_ens.push_back(evID);
  }
  return event_ID_ens;
}


void observables::output_meanpt_v2sq_correlation(){

 // create an ensemble
 std::vector<int> event_ID_ens;
 double meanpt_v2sq_corr_of_one_ens ;
 double sumx = 0 ; 
 double sumx2 = 0 ; 
  
 for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    //std::cout << ii << "  " << event_ID_ens[ii] << std::endl ; 
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_meanpt_v2sq_correlation(event_ID_ens,meanpt_v2sq_corr_of_one_ens);
    sumx += meanpt_v2sq_corr_of_one_ens ; 
    sumx2 += pow(meanpt_v2sq_corr_of_one_ens,2);  
 }
 
 // calculate mean and std. dev.
 double avg_val = sumx / rmof->get_total_music_events() ; 
 double error = sqrt( sumx2 / rmof->get_total_music_events() - pow(avg_val,2) ) ; 

 // print the output
 std::cout << "<pt>-<v_2^2> correlation = " << avg_val << "  +- " << error << std::endl ; 

}


void observables::calculate_meanpt_v2sq_correlation(std::vector<int> event_ID_ens, double& meanpt_v2sq_corr_of_one_ens){
  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ; 
  double sumv2v2star = 0. ; 
  double sumv2v2starsq = 0. ; 
  double sumptv2v2star = 0 ; 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt();
    double v2v2star = pow(ev->get_integrated_vn(2,0),2) + pow(ev->get_integrated_vn(2,1),2) ;  
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    sumv2v2star += v2v2star ;  
    sumv2v2starsq += pow(v2v2star,2) ; 
    sumptv2v2star += ( pt * v2v2star ) ;  
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  sumv2v2star /= event_ID_ens.size() ; 
  sumv2v2starsq /= event_ID_ens.size() ; 
  sumptv2v2star /= event_ID_ens.size() ; 
  double num = sumptv2v2star - sumpt * sumv2v2star ; 
  double den =  sqrt( ( sumv2v2starsq - sumv2v2star * sumv2v2star ) * ( sumptsq - sumpt * sumpt ) ) ; 
  meanpt_v2sq_corr_of_one_ens = num / den ; 
}





void observables::output_meanpt_v3sq_correlation(){

 // create an ensemble
 std::vector<int> event_ID_ens;
 double meanpt_v3sq_corr_of_one_ens ;
 double sumx = 0 ; 
 double sumx2 = 0 ; 
  
 for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    //std::cout << ii << "  " << event_ID_ens[ii] << std::endl ; 
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_meanpt_v3sq_correlation(event_ID_ens,meanpt_v3sq_corr_of_one_ens);
    sumx += meanpt_v3sq_corr_of_one_ens ; 
    sumx2 += pow(meanpt_v3sq_corr_of_one_ens,2);  
 }
 
 // calculate mean and std. dev.
 double avg_val = sumx / rmof->get_total_music_events() ; 
 double error = sqrt( sumx2 / rmof->get_total_music_events() - pow(avg_val,2) ) ; 

 // print the output
 std::cout << "<pt>-<v_3^2> correlation = " << avg_val << "  +- " << error << std::endl ; 

}


void observables::calculate_meanpt_v3sq_correlation(std::vector<int> event_ID_ens, double& meanpt_v3sq_corr_of_one_ens){
  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ; 
  double sumv3v3star = 0. ; 
  double sumv3v3starsq = 0. ; 
  double sumptv3v3star = 0 ; 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt();
    double v3v3star = pow(ev->get_integrated_vn(3,0),2) + pow(ev->get_integrated_vn(3,1),2) ;  
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    sumv3v3star += v3v3star ;  
    sumv3v3starsq += pow(v3v3star,2) ; 
    sumptv3v3star += ( pt * v3v3star ) ;  
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  sumv3v3star /= event_ID_ens.size() ; 
  sumv3v3starsq /= event_ID_ens.size() ; 
  sumptv3v3star /= event_ID_ens.size() ; 
  double num = sumptv3v3star - sumpt * sumv3v3star ; 
  double den =  sqrt( ( sumv3v3starsq - sumv3v3star * sumv3v3star ) * ( sumptsq - sumpt * sumpt ) ) ; 
  meanpt_v3sq_corr_of_one_ens = num / den ; 
}


