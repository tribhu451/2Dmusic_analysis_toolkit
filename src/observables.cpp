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

void observables::output_meanpt_vnsq_correlation(int n){

 // create an ensemble
 std::vector<int> event_ID_ens;
 
 double meanpt_vnsq_corr_of_one_ens ; // x
 double Cov_Mpt_vnvnstar; // y            
 double sumx = 0 ; 
 double sumx2 = 0 ; 
 double sumy = 0 ; 
 double sumy2 = 0 ; 
 
 double Mpt;                           
 double M_ptsq;                      
 double M_vnvnstar;                   
 double M_vnvnstarsq;                 

 // <[pt]>
 double sumA = 0. ; 
 double sumA2 = 0. ; 
 
 // (delta [pt])^2 = <[pt]^2> - <[pt]>^2
 double sumB = 0. ; 
 double sumB2 = 0. ; 
 
 // <vn * vnstar>
 double sumC = 0. ; 
 double sumC2 = 0. ; 

 // (delta vn*vnstar)^2
 double sumD = 0. ;                 
 double sumD2 = 0. ; 

 for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_meanpt_vnsq_correlation(n,event_ID_ens, Mpt, M_ptsq, M_vnvnstar, M_vnvnstarsq, Cov_Mpt_vnvnstar, meanpt_vnsq_corr_of_one_ens);
    
    sumx  += meanpt_vnsq_corr_of_one_ens ; 
    sumx2 += pow(meanpt_vnsq_corr_of_one_ens,2); 
    
    sumy  += Cov_Mpt_vnvnstar ; 
    sumy2 += pow(Cov_Mpt_vnvnstar,2); 

    sumA  += Mpt ;
    sumA2 += Mpt * Mpt ; 

    sumB  += (M_ptsq - Mpt*Mpt) ; 
    sumB2 += ( (M_ptsq - Mpt*Mpt) * (M_ptsq - Mpt*Mpt) ) ; 

    sumC  += M_vnvnstar ; 
    sumC2 += M_vnvnstar * M_vnvnstar ; 

    sumD  += ( M_vnvnstarsq - M_vnvnstar * M_vnvnstar ) ; 
    sumD2 += ( ( M_vnvnstarsq - M_vnvnstar * M_vnvnstar ) * ( M_vnvnstarsq - M_vnvnstar * M_vnvnstar ) ) ; 
 }
 
  std::ofstream mFile;
  std::stringstream output_filename;

  // calculate mean and std. dev. of rho([pt], vn*vnstar) and Cov([pt], vn*vnstar)
  double Mean_cov = sumy / rmof->get_total_music_events() ; 
  double Erro_cov = sqrt( sumy2 / rmof->get_total_music_events() - pow(Mean_cov,2) ) ; 
  double Mean_rho = sumx / rmof->get_total_music_events() ; 
  double Erro_rho = sqrt( sumx2 / rmof->get_total_music_events() - pow(Mean_rho,2) ) ; 
  // write to file
  output_filename.str("");
  output_filename << "results/Meanpt_v" << n << "v" << n << "star_correlation";
  output_filename << "_pt_";
  output_filename << ptmin << "_" << ptmax ;
  if(yflag==1){
   output_filename << "_y_" ;
  }
  else{
   output_filename << "_eta_" ;
  }
  output_filename << rapmin << "_" << rapmax ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "#Cov([pt],vn*vnstar)   error   rho([pt], vn*vnstar)   error" << std::endl ;
  mFile << Mean_cov << "   " << Erro_cov << "   " << Mean_rho << "   " << Erro_rho << std::endl ; 
  mFile.close();
  
  if(n==2){ // calculate once while calculating for n=2
    // calculate mean and std. dev. of <[pt]> and <(delta[pt])^2>
    double Mean_Mpt = sumA / rmof->get_total_music_events() ; 
    double Erro_Mpt = sqrt( sumA2 / rmof->get_total_music_events() - pow(Mean_Mpt,2) ) ; 
    double Mean_Dpt = sumB / rmof->get_total_music_events() ; 
    double Erro_Dpt = sqrt( sumB2 / rmof->get_total_music_events() - pow(Mean_Dpt,2) ) ; 
    // write to file
    output_filename.str("");
    output_filename << "results/Meanpt_and_DeltaMpt";
    output_filename << "_pt_";
    output_filename << ptmin << "_" << ptmax ;
    if(yflag==1){
     output_filename << "_y_" ;
    }
    else{
     output_filename << "_eta_" ;
    }
    output_filename << rapmin << "_" << rapmax ;
    output_filename << ".dat";
    mFile.open(output_filename.str().c_str(), std::ios::out );
    mFile << "#<[pt]>   error   <(delta [pt] )^2>   error" << std::endl ;
    mFile << Mean_Mpt << "   " << Erro_Mpt << "   " << Mean_Dpt << "   " << Erro_Dpt << std::endl ; 
    mFile.close();
  }
  
  // calculate mean and std. dev. of <vn*vnstar> and (delta vn*vnstar)^2
  double Mean_vnvnstar = sumC / rmof->get_total_music_events() ; 
  double Erro_vnvnstar = sqrt( sumC2 / rmof->get_total_music_events() - pow(Mean_vnvnstar,2) ) ; 
  double Mean_deltavnvnstar = sumD / rmof->get_total_music_events() ; 
  double Erro_deltavnvnstar = sqrt( sumD2 / rmof->get_total_music_events() - pow(Mean_deltavnvnstar,2) ) ; 
  // write to file
  output_filename.str("");
  output_filename << "results/Meanv" << n << "v" << n << "star_and_Deltav" << n << "v" << n << "star";
  output_filename << "_pt_";
  output_filename << ptmin << "_" << ptmax ;
  if(yflag==1){
   output_filename << "_y_" ;
  }
  else{
   output_filename << "_eta_" ;
  }
  output_filename << rapmin << "_" << rapmax ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "#<vn*vnstar>   error   <(delta vn*vnstar)^2>   error" << std::endl ;
  mFile << Mean_vnvnstar << "   " << Erro_vnvnstar << "   " 
  << Mean_deltavnvnstar << "   " << Erro_deltavnvnstar << std::endl ; 
  mFile.close();
}


void observables::calculate_meanpt_vnsq_correlation(int n, std::vector<int> event_ID_ens, 
     double&  Mpt, double&  M_ptsq,  double&  M_vnvnstar,  double&  M_vnvnstarsq, 
     double& Cov_Mpt_vnvnstar,   double& meanpt_vnsq_corr_of_one_ens){
  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ; 
  double sumvnvnstar = 0. ; 
  double sumvnvnstarsq = 0. ; 
  double sumptvnvnstar = 0 ; 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt();
    double vnvnstar = pow(ev->get_integrated_vn(n,0),2) + pow(ev->get_integrated_vn(n,1),2) ;  
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    sumvnvnstar += vnvnstar ;  
    sumvnvnstarsq += pow(vnvnstar,2) ; 
    sumptvnvnstar += ( pt * vnvnstar ) ;  
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  sumvnvnstar /= event_ID_ens.size() ; 
  sumvnvnstarsq /= event_ID_ens.size() ; 
  sumptvnvnstar /= event_ID_ens.size() ; 
  double num = sumptvnvnstar - sumpt * sumvnvnstar ; 
  double den =  sqrt( ( sumvnvnstarsq - sumvnvnstar * sumvnvnstar ) * ( sumptsq - sumpt * sumpt ) ) ;
  
  Mpt =  sumpt ;                             // <p_t>
  M_ptsq = sumptsq ;                         // <p_t^2>
  M_vnvnstar = sumvnvnstar ;                 // <(v_n v_n^*)>
  M_vnvnstarsq = sumvnvnstarsq ;             // <(v_n v_n^*)^2>
  Cov_Mpt_vnvnstar = num ;                   // <p_t (v_n v_n^*)> - <p_t> <(v_n v_n^*)>
  meanpt_vnsq_corr_of_one_ens = num / den ; 
}





void observables::output_pt_diff_meanpt_vnvnpt_correlation(int n){
 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> rho_meanpt_vnvnpt;
 std::vector<double> cov_meanpt_vnvnpt;
 std::vector<double> vnvnptstar;
 std::vector<double> vnvnptstarsq;

 const int ptbins = rmof->get_music_pit_bins(); 
 double sumx[ptbins];
 double sumx2[ptbins];
 double sumy[ptbins];
 double sumy2[ptbins];
 double sumA[ptbins];
 double sumA2[ptbins];
 double sumB[ptbins];
 double sumB2[ptbins];
 
 for(int ii=0; ii<ptbins; ii++){
   sumx[ii] = 0. ; 
   sumx2[ii] = 0. ; 
   sumy[ii] = 0. ; 
   sumy2[ii] = 0. ; 
   sumA[ii] = 0. ; 
   sumA2[ii] = 0. ;
   sumB[ii] = 0. ; 
   sumB2[ii] = 0. ;  
   rho_meanpt_vnvnpt.push_back(0.);
   cov_meanpt_vnvnpt.push_back(0.);
   vnvnptstar.push_back(0.);
   vnvnptstarsq.push_back(0.);
 }
  
 for(int ii=0; ii<rmof->get_total_music_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_pt_diff_meanpt_vnvnpt_correlation(n, event_ID_ens, vnvnptstar, vnvnptstarsq, cov_meanpt_vnvnpt, rho_meanpt_vnvnpt);
    for(int jj=0; jj<ptbins; jj++){
      sumx[jj]  += rho_meanpt_vnvnpt[jj] ;
      sumx2[jj] += pow(rho_meanpt_vnvnpt[jj],2) ; 
      sumy[jj]  += cov_meanpt_vnvnpt[jj] ;
      sumy2[jj] += pow(cov_meanpt_vnvnpt[jj],2) ; 
      sumA[jj]  += vnvnptstar[jj] ;
      sumA2[jj] += pow(vnvnptstar[jj],2) ; 
      sumB[jj]  += ( vnvnptstarsq[jj] - pow(vnvnptstar[jj],2) )  ;
      sumB2[jj] += ( vnvnptstarsq[jj] - pow(vnvnptstar[jj],2) ) * ( vnvnptstarsq[jj] - pow(vnvnptstar[jj],2) ) ; 
    }
 }
 
  std::ofstream mFile;
  std::stringstream output_filename;

  // calculate mean and std. dev. of cov( [pt],vn vn*(pt) ) &  rho( [pt],vn vn*(pt) )  
  // and print
  output_filename.str("");
  output_filename << "results/Meanpt_v" << n << "v" << n << "ptstar_correlation" ;
  output_filename << "_pt_";
  output_filename << ptmin << "_" << ptmax ;
  if(yflag==1){
   output_filename << "_y_" ;
  }
  else{
   output_filename << "_eta_" ;
  }
  output_filename << rapmin << "_" << rapmax ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "#pt    cov([pt],vn vn*(pt))    error    rho([pt],vn vn*(pt))    error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double cov_val = sumy[ii] / rmof->get_total_music_events() ; 
    double cov_err = sqrt( sumy2[ii] / rmof->get_total_music_events() - pow(cov_val,2) ) ; 
    double rho_val = sumx[ii] / rmof->get_total_music_events() ; 
    double rho_err = sqrt( sumx2[ii] / rmof->get_total_music_events() - pow(rho_val,2) ) ; 
    mFile << ptval << "  " << cov_val << "   " << cov_err << "   " <<  rho_val << "   " <<  rho_err << std::endl ; 
  }
  mFile.close();
 
  // calculate mean and std. dev. of ( vn vn*(pt) ) &  ( delta vn vn*(pt) )  
  // and print
  output_filename.str("");
  output_filename << "results/Meanv" << n << "v" << n << "ptstar_and_Deltav" << n << "v" << n << "ptstar" ;
  output_filename << "_pt_";
  output_filename << ptmin << "_" << ptmax ;
  if(yflag==1){
   output_filename << "_y_" ;
  }
  else{
   output_filename << "_eta_" ;
  }
  output_filename << rapmin << "_" << rapmax ;
  output_filename << ".dat";
  mFile.open(output_filename.str().c_str(), std::ios::out );
  mFile << "#pt    <vn vn*(pt)>    error    <(delta(vn vn*(pt))^2>    error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double M_val = sumA[ii] / rmof->get_total_music_events() ; 
    double M_err = sqrt( sumA2[ii] / rmof->get_total_music_events() - pow(M_val,2) ) ; 
    double D_val = sumB[ii] / rmof->get_total_music_events() ; 
    double D_err = sqrt( sumB2[ii] / rmof->get_total_music_events() - pow(D_val,2) ) ; 
    mFile << ptval << "  " << M_val << "   " << M_err << "   " <<  D_val << "   " <<  D_err << std::endl ; 
  }
  mFile.close();
}


void observables::calculate_pt_diff_meanpt_vnvnpt_correlation(int n, std::vector<int> event_ID_ens, 
   std::vector<double>& M_vnvnptstar, std::vector<double>& M_vnvnptstar_sq, 
   std::vector<double>& cov, std::vector<double>& obs){
  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ;
  const int ptbins = rmof->get_music_pit_bins(); 
  double sumvnvnptstar[ptbins] ; 
  double sumvnvnptstarsq[ptbins] ; 
  double sumptvnvnptstar[ptbins] ;
  for(int ii=0; ii<ptbins; ii++){
   sumvnvnptstar[ii] = 0. ; 
   sumvnvnptstarsq[ii] = 0. ;  
   sumptvnvnptstar[ii] = 0. ; 
  } 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt();
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    for(int jj=0; jj<ptbins; jj++){
      double vnvnptstar = ev->get_integrated_vn(n,0) * ev->get_pt_differential_vn(n,0,jj) 
       +  ev->get_integrated_vn(n,1) * ev->get_pt_differential_vn(n,1,jj) ; 
      sumvnvnptstar[jj] += vnvnptstar ; 
      sumvnvnptstarsq[jj] += pow(vnvnptstar,2) ; 
      sumptvnvnptstar[jj] += (pt * vnvnptstar) ; 
    }
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  for(int jj=0; jj<ptbins; jj++){
    sumvnvnptstar[jj] /= event_ID_ens.size() ; 
    sumvnvnptstarsq[jj] /= event_ID_ens.size() ; 
    sumptvnvnptstar[jj] /= event_ID_ens.size() ; 
  }
  
  for(int ii=0; ii<ptbins; ii++){
    M_vnvnptstar[ii] = sumvnvnptstar[ii] ; 
    M_vnvnptstar_sq[ii] = sumvnvnptstarsq[ii] ; 
    double num = sumptvnvnptstar[ii] - sumpt * sumvnvnptstar[ii] ; 
    cov[ii] = num ; 
    double den =  sqrt( ( sumvnvnptstarsq[ii] - sumvnvnptstar[ii] * sumvnvnptstar[ii] ) * ( sumptsq - sumpt * sumpt ) ) ; 
    obs[ii] = num / den ; 
  }

}











