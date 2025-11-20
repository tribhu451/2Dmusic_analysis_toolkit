#include "observables.h"
#include <cmath>

observables::observables(read_music_output_files* armof, int apid, int ayflag, double arapmin, 
   double arapmax, double aptmin, double aptmax) : rmof(armof), pid(apid), yflag(ayflag), 
     rapmin(arapmin), rapmax(arapmax), ptmin(aptmin), ptmax(aptmax){
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


void observables::output_pt_diff_meanpt_v2v2pt_correlation(){

 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> meanpt_v2v2pt_corr_of_one_ens;

 const int ptbins = rmof->get_music_pit_bins(); 
 double sumx[ptbins];
 double sumx2[ptbins];
 for(int ii=0; ii<ptbins; ii++){
   sumx[ii] = 0. ; sumx2[ii] = 0. ; 
   meanpt_v2v2pt_corr_of_one_ens.push_back(0.);
 }
  
 for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_pt_diff_meanpt_v2v2pt_correlation(event_ID_ens,meanpt_v2v2pt_corr_of_one_ens);
    for(int jj=0; jj<ptbins; jj++){
      sumx[jj] += meanpt_v2v2pt_corr_of_one_ens[jj] ;
      sumx2[jj] += pow(meanpt_v2v2pt_corr_of_one_ens[jj],2) ; 
    }
 }
 
 
  // calculate mean and std. dev. and print
  std::ofstream mFile;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "results/mpt_v2v2pt_correlation_with_pt.dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  std::cout << "=======================================" << std::endl ;
  std::cout << "calculating pt differential <pt>-v2v2pt correlation ... " << std::endl ; 
  mFile     << "#pt      correlation      error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double avg_val = sumx[ii] / rmof->get_total_music_events() ; 
    double error = sqrt( sumx2[ii] / rmof->get_total_music_events() - pow(avg_val,2) ) ; 
    //std::cout << ptval << "  " << avg_val << "   " << error << std::endl ; 
    mFile << ptval << "  " << avg_val << "   " << error << std::endl ; 
  }
  mFile.close();
 
}


void observables::calculate_pt_diff_meanpt_v2v2pt_correlation(std::vector<int> event_ID_ens, std::vector<double>& obs){

  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ;
  const int ptbins = rmof->get_music_pit_bins(); 
  double sumv2v2ptstar[ptbins] ; 
  double sumv2v2ptstarsq[ptbins] ; 
  double sumptv2v2ptstar[ptbins] ;
  for(int ii=0; ii<ptbins; ii++){
   sumv2v2ptstar[ii] = 0. ; 
   sumv2v2ptstarsq[ii] = 0. ;  
   sumptv2v2ptstar[ii] = 0. ; 
  } 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt();
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    for(int jj=0; jj<ptbins; jj++){
      double v2v2ptstar = ev->get_integrated_vn(2,0) * ev->get_pt_differential_vn(2,0,jj) 
       +  ev->get_integrated_vn(2,1) * ev->get_pt_differential_vn(2,1,jj) ; 
      sumv2v2ptstar[jj] += v2v2ptstar ; 
      sumv2v2ptstarsq[jj] += pow(v2v2ptstar,2) ; 
      sumptv2v2ptstar[jj] += (pt * v2v2ptstar) ; 
    }
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  for(int jj=0; jj<ptbins; jj++){
    sumv2v2ptstar[jj] /= event_ID_ens.size() ; 
    sumv2v2ptstarsq[jj] /= event_ID_ens.size() ; 
    sumptv2v2ptstar[jj] /= event_ID_ens.size() ; 
  }
  
  for(int ii=0; ii<ptbins; ii++){
    double num = sumptv2v2ptstar[ii] - sumpt * sumv2v2ptstar[ii] ; 
    double den =  sqrt( ( sumv2v2ptstarsq[ii] - sumv2v2ptstar[ii] * sumv2v2ptstar[ii] ) * ( sumptsq - sumpt * sumpt ) ) ; 
    obs[ii] = num / den ; 
  }

}




void observables::output_pt_diff_meanpt_v3v3pt_correlation(){

 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> meanpt_v3v3pt_corr_of_one_ens;

 const int ptbins = rmof->get_music_pit_bins(); 
 double sumx[ptbins];
 double sumx2[ptbins];
 for(int ii=0; ii<ptbins; ii++){
   sumx[ii] = 0. ; sumx2[ii] = 0. ; 
   meanpt_v3v3pt_corr_of_one_ens.push_back(0.);
 }
  
 for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_pt_diff_meanpt_v3v3pt_correlation(event_ID_ens,meanpt_v3v3pt_corr_of_one_ens);
    for(int jj=0; jj<ptbins; jj++){
      sumx[jj] += meanpt_v3v3pt_corr_of_one_ens[jj] ;
      sumx2[jj] += pow(meanpt_v3v3pt_corr_of_one_ens[jj],2) ; 
    }
 }
 
 // calculate mean and std. dev. and print
  std::ofstream mFile;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "results/mpt_v3v3pt_correlation_with_pt.dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  std::cout << "=======================================" << std::endl ;
  std::cout << "calculating pt differential <pt>-v3v3pt correlation ... " << std::endl ; 
  mFile     << "#pt      correlation      error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double avg_val = sumx[ii] / rmof->get_total_music_events() ; 
    double error = sqrt( sumx2[ii] / rmof->get_total_music_events() - pow(avg_val,2) ) ; 
    //std::cout << ptval << "  " << avg_val << "   " << error << std::endl ; 
    mFile << ptval << "  " << avg_val << "   " << error << std::endl ; 
  }
  mFile.close();
  
}


void observables::calculate_pt_diff_meanpt_v3v3pt_correlation(std::vector<int> event_ID_ens, std::vector<double>& obs){

  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ;
  const int ptbins = rmof->get_music_pit_bins(); 
  double sumv3v3ptstar[ptbins] ; 
  double sumv3v3ptstarsq[ptbins] ; 
  double sumptv3v3ptstar[ptbins] ;
  for(int ii=0; ii<ptbins; ii++){
   sumv3v3ptstar[ii] = 0. ; 
   sumv3v3ptstarsq[ii] = 0. ;  
   sumptv3v3ptstar[ii] = 0. ; 
  } 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt();
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    for(int jj=0; jj<ptbins; jj++){
      double v3v3ptstar = ev->get_integrated_vn(3,0) * ev->get_pt_differential_vn(3,0,jj) 
       +  ev->get_integrated_vn(3,1) * ev->get_pt_differential_vn(3,1,jj) ; 
      sumv3v3ptstar[jj] += v3v3ptstar ; 
      sumv3v3ptstarsq[jj] += pow(v3v3ptstar,2) ; 
      sumptv3v3ptstar[jj] += (pt * v3v3ptstar) ; 
    }
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  for(int jj=0; jj<ptbins; jj++){
    sumv3v3ptstar[jj] /= event_ID_ens.size() ; 
    sumv3v3ptstarsq[jj] /= event_ID_ens.size() ; 
    sumptv3v3ptstar[jj] /= event_ID_ens.size() ; 
  }
  
  for(int ii=0; ii<ptbins; ii++){
    double num = sumptv3v3ptstar[ii] - sumpt * sumv3v3ptstar[ii] ; 
    double den =  sqrt( ( sumv3v3ptstarsq[ii] - sumv3v3ptstar[ii] * sumv3v3ptstar[ii] ) * ( sumptsq - sumpt * sumpt ) ) ; 
    obs[ii] = num / den ; 
  }

}















