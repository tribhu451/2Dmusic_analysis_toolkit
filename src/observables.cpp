#include "observables.h"
#include <cmath>

observables::observables(read_music_output_files* armof, int ayflag, double arapmin, 
   double arapmax, double aptmin, double aptmax) : rmof(armof), yflag(ayflag), 
     rapmin(arapmin), rapmax(arapmax), ptmin(aptmin), ptmax(aptmax){
  event_arena = rmof->get_event_arena();
  rand = new random_gen();
}

std::vector<int> observables::get_an_event_ensemble(){
  std::vector<int> event_ID_ens;
  for(int ii=0; ii<rmof->get_total_events(); ii++){
    int evID =  rmof->get_total_events() * rand->rand_uniform() ;
    event_ID_ens.push_back(evID);
  }
  return event_ID_ens;
}


void observables::output_dndy_or_dndeta(){
   // species list ( charged hadron PID=0 here)
   const std::vector<int> PIDLIST =
        {0, 211,-211,321,-321,2212,-2212};
  double net_mult=0;
  std::ofstream mFile;
  std::stringstream output_filename;

  for(int PID:PIDLIST){
    net_mult = 0. ;
    for(int ii=0; ii<rmof->get_total_events(); ii++){
      int eventID = ii ; 
      event* ev = rmof->get_event(eventID) ; 
      double mult = ev->get_integrated_vn(PID,0,0); // first index is for PID, second for harmonics , third for real/imaginary
      net_mult += mult ; 
    }
    net_mult /= rmof->get_total_events(); // N
    net_mult /= (rapmax-rapmin);   // dN/deta  or dN/dy

    output_filename.str("");
    output_filename << "results/";
    if(yflag==1){
      output_filename << "dndy_";}
    else{
      output_filename << "dndeta_";}
    if(PID==0){
      output_filename << "charged_hadrons";}
    else{
      output_filename << PID;}
    output_filename << "_pt_";
    output_filename << ptmin << "_" << ptmax ;
    if(yflag==1){
     output_filename << "_y_" ; }
    else{
     output_filename << "_eta_" ;}
    output_filename << rapmin << "_" << rapmax ;
    output_filename << ".dat";
    mFile.open(output_filename.str().c_str(), std::ios::out );
    if(yflag==1){
      mFile << "#dn/dy" << std::endl ;}
    else{
      mFile << "#dn/deta" << std::endl ;
    }
    mFile <<  net_mult  << std::endl ; 
    mFile.close();
  }
  
}




void observables::output_meanpt_vnsq_correlation_charged_hadrons(int n){

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

 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_meanpt_vnsq_correlation_charged_hadrons(n,event_ID_ens, Mpt, M_ptsq, M_vnvnstar, M_vnvnstarsq, Cov_Mpt_vnvnstar, meanpt_vnsq_corr_of_one_ens);
    
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
  double Mean_cov = sumy / rmof->get_total_events() ; 
  double Erro_cov = sqrt( sumy2 / rmof->get_total_events() - pow(Mean_cov,2) ) ; 
  double Mean_rho = sumx / rmof->get_total_events() ; 
  double Erro_rho = sqrt( sumx2 / rmof->get_total_events() - pow(Mean_rho,2) ) ; 
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
    double Mean_Mpt = sumA / rmof->get_total_events() ; 
    double Erro_Mpt = sqrt( sumA2 / rmof->get_total_events() - pow(Mean_Mpt,2) ) ; 
    double Mean_Dpt = sumB / rmof->get_total_events() ; 
    double Erro_Dpt = sqrt( sumB2 / rmof->get_total_events() - pow(Mean_Dpt,2) ) ; 
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
  double Mean_vnvnstar = sumC / rmof->get_total_events() ; 
  double Erro_vnvnstar = sqrt( sumC2 / rmof->get_total_events() - pow(Mean_vnvnstar,2) ) ; 
  double Mean_deltavnvnstar = sumD / rmof->get_total_events() ; 
  double Erro_deltavnvnstar = sqrt( sumD2 / rmof->get_total_events() - pow(Mean_deltavnvnstar,2) ) ; 
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


void observables::calculate_meanpt_vnsq_correlation_charged_hadrons(int n, std::vector<int> event_ID_ens, 
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
    double pt = ev->get_mean_pt(0);
    double vnvnstar = pow(ev->get_integrated_vn(0,n,0),2) + pow(ev->get_integrated_vn(0,n,1),2) ;  
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




// prints rho([pT], v_n v_n(pT)* )  NB : The first vn is pT integrated and the other one is pT differential
void observables::output_pt_diff_meanpt_vnvnpt_correlation_charged_hadrons(int n){
 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> rho_meanpt_vnvnpt;
 std::vector<double> cov_meanpt_vnvnpt;
 std::vector<double> vnvnptstar;
 std::vector<double> vnvnptstarsq;

 const int ptbins = rmof->get_Nptbins(); 
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
  
 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_pt_diff_meanpt_vnvnpt_correlation_charged_hadrons(n, event_ID_ens, vnvnptstar, vnvnptstarsq, cov_meanpt_vnvnpt, rho_meanpt_vnvnpt);
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
    double cov_val = sumy[ii] / rmof->get_total_events() ; 
    double cov_err = sqrt( sumy2[ii] / rmof->get_total_events() - pow(cov_val,2) ) ; 
    double rho_val = sumx[ii] / rmof->get_total_events() ; 
    double rho_err = sqrt( sumx2[ii] / rmof->get_total_events() - pow(rho_val,2) ) ; 
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
    double M_val = sumA[ii] / rmof->get_total_events() ; 
    double M_err = sqrt( sumA2[ii] / rmof->get_total_events() - pow(M_val,2) ) ; 
    double D_val = sumB[ii] / rmof->get_total_events() ; 
    double D_err = sqrt( sumB2[ii] / rmof->get_total_events() - pow(D_val,2) ) ; 
    mFile << ptval << "  " << M_val << "   " << M_err << "   " <<  D_val << "   " <<  D_err << std::endl ; 
  }
  mFile.close();
}


// calculates rho([pT], v_n v_n(pT)* )  NB : The first vn is pT integrated and the other one is pT differential
void observables::calculate_pt_diff_meanpt_vnvnpt_correlation_charged_hadrons(int n, std::vector<int> event_ID_ens, 
   std::vector<double>& M_vnvnptstar, std::vector<double>& M_vnvnptstar_sq, 
   std::vector<double>& cov, std::vector<double>& obs){
  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ;
  const int ptbins = rmof->get_Nptbins(); 
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
    double pt = ev->get_mean_pt(0);
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    for(int jj=0; jj<ptbins; jj++){
      double vnvnptstar = ev->get_integrated_vn(0,n,0) * ev->get_pt_differential_vn(0,n,0,jj) 
       +  ev->get_integrated_vn(0,n,1) * ev->get_pt_differential_vn(0,n,1,jj) ; 
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




// prints rho([pT], v_n(pT) v_n(pT)* )  NB : double pT differential in vn. where in the previous one there is only one vn is pT differential.
void observables::output_pt_diff_meanpt_vnptvnpt_correlation_charged_hadrons(int n){
 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> rho_meanpt_vnptvnpt;
 std::vector<double> cov_meanpt_vnptvnpt;
 std::vector<double> vnptvnptstar;
 std::vector<double> vnptvnptstarsq;

 const int ptbins = rmof->get_Nptbins(); 
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
   rho_meanpt_vnptvnpt.push_back(0.);
   cov_meanpt_vnptvnpt.push_back(0.);
   vnptvnptstar.push_back(0.);
   vnptvnptstarsq.push_back(0.);
 }
  
 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_pt_diff_meanpt_vnptvnpt_correlation_charged_hadrons(n, event_ID_ens, vnptvnptstar, vnptvnptstarsq, cov_meanpt_vnptvnpt, rho_meanpt_vnptvnpt);
    for(int jj=0; jj<ptbins; jj++){
      sumx[jj]  += rho_meanpt_vnptvnpt[jj] ;
      sumx2[jj] += pow(rho_meanpt_vnptvnpt[jj],2) ; 
      sumy[jj]  += cov_meanpt_vnptvnpt[jj] ;
      sumy2[jj] += pow(cov_meanpt_vnptvnpt[jj],2) ; 
      sumA[jj]  += vnptvnptstar[jj] ;
      sumA2[jj] += pow(vnptvnptstar[jj],2) ; 
      sumB[jj]  += ( vnptvnptstarsq[jj] - pow(vnptvnptstar[jj],2) )  ;
      sumB2[jj] += ( vnptvnptstarsq[jj] - pow(vnptvnptstar[jj],2) ) * ( vnptvnptstarsq[jj] - pow(vnptvnptstar[jj],2) ) ; 
    }
 }
 
  std::ofstream mFile;
  std::stringstream output_filename;

  // calculate mean and std. dev. of cov( [pt],vn(pt) vn*(pt) ) &  rho( [pt],vn(pt) vn*(pt) )  
  // and print
  output_filename.str("");
  output_filename << "results/Meanpt_v" << n << "ptv" << n << "ptstar_correlation" ;
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
  mFile << "#pt    cov([pt],vn(pt) vn*(pt))    error    rho([pt],vn(pt) vn*(pt))    error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double cov_val = sumy[ii] / rmof->get_total_events() ; 
    double cov_err = sqrt( sumy2[ii] / rmof->get_total_events() - pow(cov_val,2) ) ; 
    double rho_val = sumx[ii] / rmof->get_total_events() ; 
    double rho_err = sqrt( sumx2[ii] / rmof->get_total_events() - pow(rho_val,2) ) ; 
    mFile << ptval << "  " << cov_val << "   " << cov_err << "   " <<  rho_val << "   " <<  rho_err << std::endl ; 
  }
  mFile.close();
 
  // calculate mean and std. dev. of ( vn(pt) vn*(pt) ) &  ( delta vn(pt) vn*(pt) )  
  // and print
  output_filename.str("");
  output_filename << "results/Meanv" << n << "ptv" << n << "ptstar_and_Deltav" << n << "ptv" << n << "ptstar" ;
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
  mFile << "#pt    <vn(pt) vn*(pt)>    error    <(delta(vn(pt) vn*(pt))^2>    error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double M_val = sumA[ii] / rmof->get_total_events() ; 
    double M_err = sqrt( sumA2[ii] / rmof->get_total_events() - pow(M_val,2) ) ; 
    double D_val = sumB[ii] / rmof->get_total_events() ; 
    double D_err = sqrt( sumB2[ii] / rmof->get_total_events() - pow(D_val,2) ) ; 
    mFile << ptval << "  " << M_val << "   " << M_err << "   " <<  D_val << "   " <<  D_err << std::endl ; 
  }
  mFile.close();
}




// calculates rho([pT], v_n(pT) v_n(pT)* )  NB : double pT differential in vn. where in the previous one there is only one vn is pT differential.
void observables::calculate_pt_diff_meanpt_vnptvnpt_correlation_charged_hadrons(int n, std::vector<int> event_ID_ens, 
   std::vector<double>& M_vnptvnptstar, std::vector<double>& M_vnptvnptstar_sq, 
   std::vector<double>& cov, std::vector<double>& obs){
  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ;
  const int ptbins = rmof->get_Nptbins(); 
  double sumvnptvnptstar[ptbins] ; 
  double sumvnptvnptstarsq[ptbins] ; 
  double sumptvnptvnptstar[ptbins] ;
  for(int ii=0; ii<ptbins; ii++){
   sumvnptvnptstar[ii] = 0. ; 
   sumvnptvnptstarsq[ii] = 0. ;  
   sumptvnptvnptstar[ii] = 0. ; 
  } 
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt(0);
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    for(int jj=0; jj<ptbins; jj++){
      double vnptvnptstar = ev->get_pt_differential_vn(0,n,0,jj)  * ev->get_pt_differential_vn(0,n,0,jj) 
       +   ev->get_pt_differential_vn(0,n,1,jj) * ev->get_pt_differential_vn(0,n,1,jj) ; 
      sumvnptvnptstar[jj] += vnptvnptstar ; 
      sumvnptvnptstarsq[jj] += pow(vnptvnptstar,2) ; 
      sumptvnptvnptstar[jj] += (pt * vnptvnptstar) ; 
    }
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  for(int jj=0; jj<ptbins; jj++){
    sumvnptvnptstar[jj] /= event_ID_ens.size() ; 
    sumvnptvnptstarsq[jj] /= event_ID_ens.size() ; 
    sumptvnptvnptstar[jj] /= event_ID_ens.size() ; 
  }
  
  for(int ii=0; ii<ptbins; ii++){
    M_vnptvnptstar[ii] = sumvnptvnptstar[ii] ; 
    M_vnptvnptstar_sq[ii] = sumvnptvnptstarsq[ii] ; 
    double num = sumptvnptvnptstar[ii] - sumpt * sumvnptvnptstar[ii] ; 
    cov[ii] = num ; 
    double den =  sqrt( ( sumvnptvnptstarsq[ii] - sumvnptvnptstar[ii] * sumvnptvnptstar[ii] ) * ( sumptsq - sumpt * sumpt ) ) ; 
    obs[ii] = num / den ; 
  }

}



void observables::output_pt_diff_multiparticle_vn_charged_hadrons(int n){
 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> vn_2_pt;
 std::vector<double> vn_4_pt;
 
 double vn_2part_sq;
 double vn_4part_fr;

 const int ptbins = rmof->get_Nptbins(); 
 double sumx[ptbins];
 double sumx2[ptbins];
 double sumy[ptbins];
 double sumy2[ptbins];

 double sumw1 = 0 ; 
 double sumw1_sq = 0 ; 
 double sumw2 = 0 ; 
 double sumw2_sq = 0 ; 
 
 for(int ii=0; ii<ptbins; ii++){
   sumx[ii] = 0. ; 
   sumx2[ii] = 0. ; 
   sumy[ii] = 0. ; 
   sumy2[ii] = 0. ; 
   vn_2_pt.push_back(0.);
   vn_4_pt.push_back(0.);
 }
  
   //for(int ii=0; ii<rmof->get_total_events(); ii++){
   int iEns = 0 ;
   int itry = 0 ;
   do{
     event_ID_ens = get_an_event_ensemble();
     itry++ ; 
     // calculate correlation of a given ensemble
     calculate_pt_diff_multiparticle_vn_charged_hadrons(n, event_ID_ens, vn_2part_sq, vn_4part_fr,  vn_2_pt, vn_4_pt );
     if(vn_2part_sq < 0 || vn_4part_fr < 0 ){
       continue ; 
     }else{
       iEns++ ; 
       sumw1 += sqrt(vn_2part_sq) ;  //  vn{2}
       sumw1_sq += vn_2part_sq  ; 
       sumw2 += pow(vn_4part_fr,0.25) ;  // vn{4}
       sumw2_sq += pow(vn_4part_fr,0.5) ; 
       for(int jj=0; jj<ptbins; jj++){
        sumx[jj]  += vn_2_pt[jj] ;
        sumx2[jj] += pow(vn_2_pt[jj],2) ; 
        sumy[jj]  += vn_4_pt[jj] ;
        sumy2[jj] += pow(vn_4_pt[jj],2) ; 
       }
     }
   }
   while(iEns<rmof->get_total_events());
   std::cout << "Number of ensemble tries to calculate v"<<n<<"{2} and {4} = " << itry << std::endl;
 
  std::ofstream mFile;
  std::stringstream output_filename;

  // calculate mean and std. dev. 
  // and print
  output_filename.str("");
  output_filename << "results/Multiparticle_v" << n ;
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
  mFile << "#pt  int-v" << n <<"{2}   error   int-v" << n << "{4}    error" << "  v" << n<<"{2}(pt)    error    v_"<<n<<"{4}(pt)    error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);
    double int_v2_2 = sumw1 / rmof->get_total_events() ; 
    double int_v2_2_err = sqrt( sumw1_sq / rmof->get_total_events() - pow(int_v2_2,2) ) ;    
    double int_v2_4 = sumw2 / rmof->get_total_events() ; 
    double int_v2_4_err = sqrt( sumw2_sq / rmof->get_total_events() - pow(int_v2_4,2) ) ;    
    double vn_2_val = sumx[ii] / rmof->get_total_events() ; 
    double vn_2_err = sqrt( sumx2[ii] / rmof->get_total_events() - pow(vn_2_val,2) ) ; 
    double vn_4_val = sumy[ii] / rmof->get_total_events() ; 
    double vn_4_err = sqrt( sumy2[ii] / rmof->get_total_events() - pow(vn_4_val,2) ) ; 
    mFile << ptval << "  " << int_v2_2 << "  " << int_v2_2_err << "  " << int_v2_4 << "  " <<  int_v2_4_err 
    << "  " << vn_2_val << "   " << vn_2_err << "   " <<  vn_4_val << "   " <<  vn_4_err << std::endl ; 
  }
  mFile.close();

}



// calculates v_n{2}(pT) and v_n{4}(pT)
void observables::calculate_pt_diff_multiparticle_vn_charged_hadrons(int n, std::vector<int> event_ID_ens, 
   double& vn_sq, double& vn_fr, std::vector<double>& vn_2, std::vector<double>& vn_4){
  // calculate the observable for one ensemble //
  const int ptbins = rmof->get_Nptbins(); 
  double sum1[ptbins] ; 
  double sum2 = 0. ; 
  double sum3[ptbins] ;
  double sum4 = 0. ;
  for(int ii=0; ii<ptbins; ii++){
   sum1[ii] = 0. ; 
   sum3[ii] = 0. ; 
  } 
  double Cn, Sn, Cn_pt, Sn_pt, Cn_sq, Sn_sq ; 
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    Cn    = ev->get_integrated_vn(0,n,0) ; 
    Sn    = ev->get_integrated_vn(0,n,1) ;
    Cn_sq = Cn * Cn ; 
    Sn_sq = Sn * Sn ; 
    sum2 += ( Cn_sq + Sn_sq ) ; // < v2 v2*>
    sum4 += ( ( Cn_sq + Sn_sq ) * ( Cn_sq + Sn_sq ) ) ; // < v2 v2* v2 v2* >
    for(int jj=0; jj<ptbins; jj++){
      Cn_pt = ev->get_pt_differential_vn(0,n,0,jj) ;
      Sn_pt = ev->get_pt_differential_vn(0,n,1,jj) ;
      sum1[jj] += ( Cn_sq + Sn_sq ) * ( Cn * Cn_pt + Sn * Sn_pt ) ; // < v2 v2* v2 v2*(pt) >
      sum3[jj] += ( Cn * Cn_pt + Sn * Sn_pt ) ; // < v2 v2*(pt) >
    }
  }
  
  sum2 /= event_ID_ens.size() ; 
  sum4 /= event_ID_ens.size() ; 

  for(int jj=0; jj<ptbins; jj++){
    sum1[jj] /= event_ID_ens.size() ; 
    sum3[jj] /= event_ID_ens.size() ; 
  }
  
  double den = 2. * sum2 * sum2 - sum4   ;

  vn_sq = sum2 ; // < v2 v2*>
  vn_fr = den ;  // 2 < v2 v2*> < v2 v2*> - < v2 v2* v2 v2* > 
    
  for(int ii=0; ii<ptbins; ii++){
    if(vn_sq < 0 || vn_fr < 0 ){
     vn_2[ii]=0.;
     vn_4[ii]=0.;
    }
    else{
      vn_2[ii] = sum3[ii] / sqrt(sum2) ;
      vn_4[ii] = ( 2 * sum2 * sum3[ii] - sum1[ii]) / pow(den,3./4.) ; 
    }
  }

}



// In method 2,  the numerator and denominator of vn{2} nad vn{4}
// are calculated independently and then their error is calculated using 
// quadrature sum of errors of both numerator and denominator.
void observables::output_pt_diff_multiparticle_vn_method2_charged_hadrons(int n){
 // create an ensemble
 std::vector<int> event_ID_ens;
 std::vector<double> vn_2_num_pt;
 std::vector<double> vn_4_num_pt;
 
 double vn_2part_sq;
 double vn_4part_fr;

 const int ptbins = rmof->get_Nptbins(); 
 double sumx[ptbins];
 double sumx2[ptbins];
 double sumy[ptbins];
 double sumy2[ptbins];

 double sumw1 = 0 ; 
 double sumw1_sq = 0 ; 
 double sumw2 = 0 ; 
 double sumw2_sq = 0 ; 
 
 for(int ii=0; ii<ptbins; ii++){
   sumx[ii] = 0. ; 
   sumx2[ii] = 0. ; 
   sumy[ii] = 0. ; 
   sumy2[ii] = 0. ; 
   vn_2_num_pt.push_back(0.);
   vn_4_num_pt.push_back(0.);
 }
  
   //for(int ii=0; ii<rmof->get_total_events(); ii++){
   int iEns = 0 ;
   do{
     event_ID_ens = get_an_event_ensemble();
     // calculate vn of a given ensemble
     calculate_pt_diff_multiparticle_vn_method2_charged_hadrons(n, event_ID_ens, vn_2part_sq, vn_4part_fr, vn_2_num_pt, vn_4_num_pt );
     iEns++ ; 
     sumw1 += vn_2part_sq ;  // < v2 v2*>
     sumw1_sq += vn_2part_sq * vn_2part_sq  ; 
     sumw2 += vn_4part_fr ;  // 2 < v2 v2*> < v2 v2*> - < v2 v2* v2 v2* > 
     sumw2_sq += vn_4part_fr * vn_4part_fr ; 
     for(int jj=0; jj<ptbins; jj++){
      sumx[jj]  += vn_2_num_pt[jj] ;
      sumx2[jj] += pow(vn_2_num_pt[jj],2) ; 
      sumy[jj]  += vn_4_num_pt[jj] ;
      sumy2[jj] += pow(vn_4_num_pt[jj],2) ;    
     }
   }
   while(iEns<rmof->get_total_events());
   
   sumw1 /=  rmof->get_total_events() ; 
   sumw2 /=  rmof->get_total_events() ; 
   sumw1_sq /=  rmof->get_total_events() ; 
   sumw2_sq /=  rmof->get_total_events() ; 
   if(sumw1<0){
     std::cout << "error in method2 calculation of vn{2}, vn{4} ..." << std::endl ; 
     std::cout << "<v2 v2*> negative ..." << std::endl ; 
     exit(-1);
   }
   if(sumw2<0){
     std::cout << "error in method2 calculation of vn{2}, vn{4} ..." << std::endl ; 
     std::cout << "( 2 < v2 v2*> < v2 v2*> - < v2 v2* v2 v2* > )  negative ..." << std::endl ; 
     exit(-1);
   }
   
   // now calculate the integrated vn{2},vn{4} and it's error.
   double temp1, temp2;
   double inti_vn_2, inti_vn_2_err, inti_vn_4, inti_vn_4_err;
   //vn{2}
   inti_vn_2 = sqrt(sumw1);
   temp1 = sqrt(sumw1_sq - sumw1 * sumw1) ;
   inti_vn_2_err = temp1 * 0.50 * 1. / sqrt(sumw1);
   // vn{4}
   inti_vn_4 = pow(sumw2,0.25);
   temp1 = sqrt(sumw2_sq - sumw2 * sumw2) ;
   inti_vn_4_err = temp1 * 0.25 * pow(sumw2,-0.75);
   
   for(int ii=0; ii<ptbins; ii++){
     sumx[ii] /= rmof->get_total_events() ; 
     sumx2[ii] /= rmof->get_total_events() ; 
     sumy[ii] /= rmof->get_total_events() ; 
     sumy2[ii] /= rmof->get_total_events() ; 
   }   

   // now calculate the pt-differential vn{2},vn{4} and it's error.
   double vn_2[ptbins];
   double vn_4[ptbins];
   double vn_2_err[ptbins];
   double vn_4_err[ptbins];
   for(int ii=0; ii<ptbins; ii++){ 
     vn_2[ii] = sumx[ii] / inti_vn_2 ;
     temp1 = sqrt(sumx2[ii] - sumx[ii] * sumx[ii]);
     vn_2_err[ii] = sqrt( pow(sumx[ii],2)/pow(inti_vn_2,4) * pow(inti_vn_2_err,2) + 1. / pow(inti_vn_2,2) * pow(temp1,2) ) ;
     
     vn_4[ii] = sumy[ii] / pow(sumw2,0.75); ; // y = sumy[ii], Dy = temp1,  x = sumw2, Dx = temp2 
     temp1 = sqrt(sumy2[ii] - sumy[ii] * sumy[ii]);
     temp2 =  sqrt(sumw2_sq - sumw2 * sumw2) ;
     vn_4_err[ii] = sqrt( 9./16. * pow(sumy[ii],2) * pow(sumw2,-7./2.) * pow(temp2,2) + pow(sumw2,-3./2.)*pow(temp1,2) ) ;
   }


  std::ofstream mFile;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "results/Multiparticle_v" << n << "_method2";
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
  mFile << "#pt  int-v" << n <<"{2}   error   int-v" << n << "{4}    error" << "  v" << n<<"{2}(pt)    error    v_"<<n<<"{4}(pt)    error" << std::endl ; 
  for(int ii=0; ii<ptbins; ii++){
    double ptval = rmof->get_pt_val_of_bin(ii);  
    mFile << ptval << "  " << inti_vn_2 << "  " << inti_vn_2_err << "  " << inti_vn_4 << "  " <<  inti_vn_4_err <<
    "  " << vn_2[ii] << "  " << vn_2_err[ii] << "  " << vn_4[ii] << "  " << vn_4_err[ii] << std::endl ; 
  }
  mFile.close();

}



// calculates v_n{2}(pT) and v_n{4}(pT)
void observables::calculate_pt_diff_multiparticle_vn_method2_charged_hadrons(int n, std::vector<int> event_ID_ens, 
   double& vn_sq, double& vn_fr, std::vector<double>& vn_2_num, std::vector<double>& vn_4_num ){
  // calculate the observable for one ensemble //
  const int ptbins = rmof->get_Nptbins(); 
  double sum1[ptbins] ; 
  double sum2 = 0. ; 
  double sum3[ptbins] ;
  double sum4 = 0. ;
  for(int ii=0; ii<ptbins; ii++){
   sum1[ii] = 0. ; 
   sum3[ii] = 0. ; 
  } 
  double Cn, Sn, Cn_pt, Sn_pt, Cn_sq, Sn_sq ; 
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    Cn    = ev->get_integrated_vn(0,n,0) ; 
    Sn    = ev->get_integrated_vn(0,n,1) ;
    Cn_sq = Cn * Cn ; 
    Sn_sq = Sn * Sn ; 
    sum2 += ( Cn_sq + Sn_sq ) ; // < v2 v2*>
    sum4 += ( ( Cn_sq + Sn_sq ) * ( Cn_sq + Sn_sq ) ) ; // < v2 v2* v2 v2* >
    for(int jj=0; jj<ptbins; jj++){
      Cn_pt = ev->get_pt_differential_vn(0,n,0,jj) ;
      Sn_pt = ev->get_pt_differential_vn(0,n,1,jj) ;
      sum1[jj] += ( Cn_sq + Sn_sq ) * ( Cn * Cn_pt + Sn * Sn_pt ) ; // < v2 v2* v2 v2*(pt) >
      sum3[jj] += ( Cn * Cn_pt + Sn * Sn_pt ) ; // < v2 v2*(pt) >
    }
  }
  
  sum2 /= event_ID_ens.size() ; 
  sum4 /= event_ID_ens.size() ; 

  for(int jj=0; jj<ptbins; jj++){
    sum1[jj] /= event_ID_ens.size() ; 
    sum3[jj] /= event_ID_ens.size() ; 
  }
  
  double den = 2. * sum2 * sum2 - sum4   ;
  vn_sq = sum2 ; // < v2 v2*>
  vn_fr = den ;  // 2 < v2 v2*> < v2 v2*> - < v2 v2* v2 v2* > 
    
  for(int ii=0; ii<ptbins; ii++){
   vn_2_num[ii] = sum3[ii] ;
   vn_4_num[ii] = ( 2. * sum2 * sum3[ii] - sum1[ii] ) ; 
  }

}







void observables::output_meanpt_vnsq_higher_moments_charged_hadrons(int n){

 // create an ensemble
 std::vector<int> event_ID_ens;
 
 double sumA = 0 ; 
 double sumAsq = 0 ; 
 double sumB = 0 ; 
 double sumBsq = 0 ; 
 double sumC = 0 ; 
 double sumCsq = 0 ; 
 double sumD = 0 ; 
 double sumDsq = 0 ; 
 double sumE = 0 ; 
 double sumEsq = 0 ; 

 double sumU = 0 ; 
 double sumUsq = 0 ; 
 double sumV = 0 ; 
 double sumVsq = 0 ; 
 double sumW = 0 ; 
 double sumWsq = 0 ; 
 double sumX = 0 ; 
 double sumXsq = 0 ; 
 double sumY = 0 ; 
 double sumYsq = 0 ; 

 double c11, c21, c31, c41, c122;
 double r11, r21, r31, r41, r122;
 
 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_meanpt_vnsq_higher_moments_charged_hadrons(n,event_ID_ens, c11, r11, c21, r21, c31, r31, c41, r41, c122, r122);
    sumA  += c11 ; 
    sumAsq += pow(c11,2); 
    sumB  += c21 ; 
    sumBsq += pow(c21,2); 
    sumC  += c31 ; 
    sumCsq += pow(c31,2); 
    sumD  += c41 ; 
    sumDsq += pow(c41,2); 
    sumE  += c122 ; 
    sumEsq += pow(c122,2); 
    
    sumU  += r11 ; 
    sumUsq += pow(r11,2); 
    sumV  += r21 ; 
    sumVsq += pow(r21,2); 
    sumW  += r31 ; 
    sumWsq += pow(r31,2); 
    sumX  += r41 ; 
    sumXsq += pow(r41,2); 
    sumY  += r122 ; 
    sumYsq += pow(r122,2); 
 }
 
  
  std::ofstream mFile;
  std::stringstream output_filename;

  double val;
  double err;

  // write to file
  output_filename.str("");
  output_filename << "results/Meanpt_v" << n << "sq_higher_moments";
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
  mFile << "#Cov(1,1)   error   Cov(2,1)   error   Cov(3,1)   error    Cov(4,1)   error  Cov(1,2,2)   error " 
  << "  rho(1,1)    error    rho(2,1)    error   rho(3,1)    error  rho(4,1)  error   rho(1,2,2)   error  " 
  << std::endl ;

  val = sumA / rmof->get_total_events() ; 
  err = sqrt( sumAsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  "; 
  
  val = sumB / rmof->get_total_events() ; 
  err = sqrt( sumBsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumC / rmof->get_total_events() ; 
  err = sqrt( sumCsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumD / rmof->get_total_events() ; 
  err = sqrt( sumDsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumE / rmof->get_total_events() ; 
  err = sqrt( sumEsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 
 
  val = sumU / rmof->get_total_events() ; 
  err = sqrt( sumUsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 
  
  val = sumV / rmof->get_total_events() ; 
  err = sqrt( sumVsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumW / rmof->get_total_events() ; 
  err = sqrt( sumWsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumX / rmof->get_total_events() ; 
  err = sqrt( sumXsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumY / rmof->get_total_events() ; 
  err = sqrt( sumYsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 
   
  mFile << std::endl ; 
  mFile.close();
  
}



void observables::calculate_meanpt_vnsq_higher_moments_charged_hadrons(int n, std::vector<int> event_ID_ens, 
    double& Cov11, double& rho11, double& Cov21, double& rho21, 
    double& Cov31, double& rho31, double& Cov41, double& rho41, 
    double& Cov122, double& rho122 ){

  // calculate the observable for one ensemble //
  double sumpt = 0. ; 
  double sumptsq = 0. ; 
  double sumvnvnstar = 0. ; 
  double sumvnvnstarsq = 0. ; 

  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt(0);
    double vnvnstar = pow(ev->get_integrated_vn(0,n,0),2) + pow(ev->get_integrated_vn(0,n,1),2) ;  
    sumpt += pt ; 
    sumptsq += pow(pt,2);
    sumvnvnstar += vnvnstar ;  
    sumvnvnstarsq += pow(vnvnstar,2) ; 
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumptsq /= event_ID_ens.size() ; 
  sumvnvnstar /= event_ID_ens.size() ; 
  sumvnvnstarsq /= event_ID_ens.size() ; 

  double avg_mpt  = sumpt ; 
  double avg_vnsq = sumvnvnstar ; 
  
  double sigma_pt = sqrt( ( sumptsq - sumpt * sumpt ) ) ;
  double sigma_vnsq = sqrt( ( sumvnvnstarsq - sumvnvnstar * sumvnvnstar )  ) ;
  
  double sum_dptsq=0.; 
  double sum_dptcu=0.; 
  double sum_dpt_dvnsq=0.; 
  double sum_dptsq_dvnsq=0.; 
  double sum_dptcu_dvnsq=0.; 
  double sum_dpt4_dvnsq=0.; 
  double sum_dpt_dvnsq_dvnsq=0.;
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double pt = ev->get_mean_pt(0);
    double vnvnstar = pow(ev->get_integrated_vn(0,n,0),2) + pow(ev->get_integrated_vn(0,n,1),2) ;  
    sum_dptsq += pow(pt-avg_mpt,2.);
    sum_dptcu += pow(pt-avg_mpt,3.);
    sum_dpt_dvnsq +=  (pt-avg_mpt)*(vnvnstar - avg_vnsq);
    sum_dptsq_dvnsq +=  pow(pt-avg_mpt,2.)*(vnvnstar - avg_vnsq);
    sum_dptcu_dvnsq +=  pow(pt-avg_mpt,3.)*(vnvnstar - avg_vnsq);
    sum_dpt4_dvnsq  +=  pow(pt-avg_mpt,4.)*(vnvnstar - avg_vnsq);
    sum_dpt_dvnsq_dvnsq += (pt-avg_mpt)*(vnvnstar - avg_vnsq)*(vnvnstar - avg_vnsq);
  }
  
  sum_dptsq /= event_ID_ens.size() ; 
  sum_dptcu /= event_ID_ens.size() ; 
  sum_dpt_dvnsq /= event_ID_ens.size() ; 
  sum_dptsq_dvnsq /= event_ID_ens.size() ; 
  sum_dptcu_dvnsq /= event_ID_ens.size() ; 
  sum_dpt4_dvnsq /= event_ID_ens.size() ; 
  sum_dpt_dvnsq_dvnsq /= event_ID_ens.size() ; 

  double num ; 
  double den ;
   
  num = sum_dpt_dvnsq ; 
  den = pow(sigma_pt,1) * sigma_vnsq ;
  Cov11 = num  ; 
  rho11 = num / den ; 
  
  num = sum_dptsq_dvnsq ; 
  den = pow(sigma_pt,2) * sigma_vnsq ;
  Cov21 = num  ; 
  rho21 = num / den ; 

  num = sum_dptcu_dvnsq - 3 * sum_dptsq * sum_dpt_dvnsq ; 
  den = pow(sigma_pt,3) * sigma_vnsq ;
  Cov31 = num  ; 
  rho31 = num / den ; 

  num = sum_dpt4_dvnsq - 6 * sum_dptsq *  sum_dptsq_dvnsq - 4 * sum_dptcu * sum_dpt_dvnsq ; 
  den = pow(sigma_pt,4) * sigma_vnsq ;
  Cov41 = num  ; 
  rho41 = num / den ; 

  num = sum_dpt_dvnsq_dvnsq - 2 * avg_vnsq * sum_dpt_dvnsq ; 
  den = pow(sigma_pt,1) * pow(sigma_vnsq,2);
  Cov122 = num  ; 
  rho122 = num / den ; 

}






void observables::output_meanpt_vnsq_higher_moments_mult_fluc_corrected_charged_hadrons(int n){

 // create an ensemble
 std::vector<int> event_ID_ens;
 
 double sumA = 0 ; 
 double sumAsq = 0 ; 
 double sumB = 0 ; 
 double sumBsq = 0 ; 
 double sumC = 0 ; 
 double sumCsq = 0 ; 
 double sumD = 0 ; 
 double sumDsq = 0 ; 
 double sumE = 0 ; 
 double sumEsq = 0 ; 

 double sumU = 0 ; 
 double sumUsq = 0 ; 
 double sumV = 0 ; 
 double sumVsq = 0 ; 
 double sumW = 0 ; 
 double sumWsq = 0 ; 
 double sumX = 0 ; 
 double sumXsq = 0 ; 
 double sumY = 0 ; 
 double sumYsq = 0 ; 

 double c11, c21, c31, c41, c122;
 double r11, r21, r31, r41, r122;
 
 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_meanpt_vnsq_higher_moments_mult_fluc_corrected_charged_hadrons(n,event_ID_ens, c11, r11, c21, r21, c31, r31, c41, r41, c122, r122);
    sumA  += c11 ; 
    sumAsq += pow(c11,2); 
    sumB  += c21 ; 
    sumBsq += pow(c21,2); 
    sumC  += c31 ; 
    sumCsq += pow(c31,2); 
    sumD  += c41 ; 
    sumDsq += pow(c41,2); 
    sumE  += c122 ; 
    sumEsq += pow(c122,2); 
    
    sumU  += r11 ; 
    sumUsq += pow(r11,2); 
    sumV  += r21 ; 
    sumVsq += pow(r21,2); 
    sumW  += r31 ; 
    sumWsq += pow(r31,2); 
    sumX  += r41 ; 
    sumXsq += pow(r41,2); 
    sumY  += r122 ; 
    sumYsq += pow(r122,2); 
 }
 
  
  std::ofstream mFile;
  std::stringstream output_filename;

  double val;
  double err;

  // write to file
  output_filename.str("");
  output_filename << "results/Meanpt_v" << n << "sq_higher_moments_mult_fluc_corrected";
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
  mFile << "#Cov(1,1)   error   Cov(2,1)   error   Cov(3,1)   error    Cov(4,1)   error  Cov(1,2,2)   error " 
  << "  rho(1,1)    error    rho(2,1)    error   rho(3,1)    error  rho(4,1)  error   rho(1,2,2)   error  " 
  << std::endl ;

  val = sumA / rmof->get_total_events() ; 
  err = sqrt( sumAsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  "; 
  
  val = sumB / rmof->get_total_events() ; 
  err = sqrt( sumBsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumC / rmof->get_total_events() ; 
  err = sqrt( sumCsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumD / rmof->get_total_events() ; 
  err = sqrt( sumDsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumE / rmof->get_total_events() ; 
  err = sqrt( sumEsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 
 
  val = sumU / rmof->get_total_events() ; 
  err = sqrt( sumUsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 
  
  val = sumV / rmof->get_total_events() ; 
  err = sqrt( sumVsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumW / rmof->get_total_events() ; 
  err = sqrt( sumWsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumX / rmof->get_total_events() ; 
  err = sqrt( sumXsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 

  val = sumY / rmof->get_total_events() ; 
  err = sqrt( sumYsq / rmof->get_total_events() - pow(val,2) ) ; 
  mFile << val << "  " << err << "  " ; 
   
  mFile << std::endl ; 
  mFile.close();
  
}


void observables::calculate_meanpt_vnsq_higher_moments_mult_fluc_corrected_charged_hadrons(int n, std::vector<int> event_ID_ens, 
    double& Cov11, double& rho11, double& Cov21, double& rho21, 
    double& Cov31, double& rho31, double& Cov41, double& rho41, 
    double& Cov122, double& rho122 ){

   // frequently used variables, globally declare karagala.
   // allocation inside loop takes time.
   double nch;
   double spt;
   double vnvnstar;
   double deltapt;
   double deltavnsq;


  // first : Var(Nch), average <pt> and average <vn^2> //
  double sumpt = 0. ; 
  double sumnch = 0. ; 
  double sumnchsq = 0. ; 
  double sumvnvnstar = 0. ; 

  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    nch = ev->get_integrated_vn(0,0,0);
    spt = ev->get_mean_pt(0);
    vnvnstar = pow(ev->get_integrated_vn(0,n,0),2) + pow(ev->get_integrated_vn(0,n,1),2) ;  
    sumpt += spt ; 
    sumvnvnstar += vnvnstar ;  
    sumnch += nch ; 
    sumnchsq += pow(nch,2);
  }
  
  sumpt /= event_ID_ens.size() ; 
  sumvnvnstar /= event_ID_ens.size() ; 
  sumnch /= event_ID_ens.size() ; 
  sumnchsq /= event_ID_ens.size() ; 
  double avg_nch = sumnch ;                                      // <Nch>
  double variance_nch =  ( sumnchsq - sumnch * sumnch )  ;       // Var(Nch)
  double avg_mpt = sumpt ;                                       // <pt>
  double avg_vnsq = sumvnvnstar ;                                // <vn^2>


  // second : cov(delta pt, Nch) and cov(delta vn^2, Nch) //
  double sum_deltapt_nch = 0 ; 
  double sum_deltavnsq_nch = 0 ; 
  double sum_deltapt = 0 ; 
  double sum_deltavnsq = 0 ;
  double sum_nch = 0 ;  
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    nch = ev->get_integrated_vn(0,0,0);
    spt = ev->get_mean_pt(0);
    vnvnstar = pow(ev->get_integrated_vn(0,n,0),2) + pow(ev->get_integrated_vn(0,n,1),2) ; 
    // deltapt 
    deltapt = (spt - avg_mpt) ;
    // delta vn^2  
    deltavnsq = (vnvnstar - avg_vnsq) ;  
    
    sum_deltapt_nch   += (deltapt*nch);
    sum_deltavnsq_nch += (deltavnsq*nch);
    sum_deltapt       += (deltapt);
    sum_deltavnsq     += (deltavnsq);
    sum_nch           += (nch);
  }
  
  sum_deltapt_nch   /= event_ID_ens.size() ; 
  sum_deltavnsq_nch /= event_ID_ens.size() ; 
  sum_deltapt       /= event_ID_ens.size() ; 
  sum_deltavnsq     /= event_ID_ens.size() ; 
  sum_nch           /= event_ID_ens.size() ; 

  double cov_deltapt_nch = sum_deltapt_nch - sum_deltapt * sum_nch ;        // Cov(delta pt, Nch)
  double cov_deltavnsq_nch = sum_deltavnsq_nch - sum_deltavnsq * sum_nch ;  // Cov(delta vn^2, Nch)


  // multiplicity fluctuation corrected calculation //
  double sum_dptsq=0.; 
  double sum_dvn2sq=0.; 
  double sum_dptcu=0.; 
  double sum_dpt_dvnsq=0.; 
  double sum_dptsq_dvnsq=0.; 
  double sum_dptcu_dvnsq=0.; 
  double sum_dpt4_dvnsq=0.; 
  double sum_dpt_dvnsq_dvnsq=0.;
   
  double deltapt_tilde, delta_vnsq_tilde;
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    spt = ev->get_mean_pt(0);
    nch = ev->get_integrated_vn(0,0,0);
    vnvnstar = pow(ev->get_integrated_vn(0,n,0),2) + pow(ev->get_integrated_vn(0,n,1),2) ; 
    deltapt_tilde =  (spt-avg_mpt) - cov_deltapt_nch / variance_nch * ( nch - avg_nch ) ; 
    delta_vnsq_tilde =  (vnvnstar-avg_vnsq) - cov_deltavnsq_nch / variance_nch * ( nch - avg_nch ) ; 
    sum_dvn2sq += pow(delta_vnsq_tilde,2.);
    sum_dptsq += pow(deltapt_tilde,2.);
    sum_dptcu += pow(deltapt_tilde,3.);
    sum_dpt_dvnsq +=  (deltapt_tilde)*(delta_vnsq_tilde);
    sum_dptsq_dvnsq +=  pow(deltapt_tilde,2.)*(delta_vnsq_tilde);
    sum_dptcu_dvnsq +=  pow(deltapt_tilde,3.)*(delta_vnsq_tilde);
    sum_dpt4_dvnsq  +=  pow(deltapt_tilde,4.)*(delta_vnsq_tilde);
    sum_dpt_dvnsq_dvnsq += (deltapt_tilde)*(delta_vnsq_tilde)*(delta_vnsq_tilde);
  }
  
  sum_dvn2sq /= event_ID_ens.size() ; 
  sum_dptsq /= event_ID_ens.size() ; 
  sum_dptcu /= event_ID_ens.size() ; 
  sum_dpt_dvnsq /= event_ID_ens.size() ; 
  sum_dptsq_dvnsq /= event_ID_ens.size() ; 
  sum_dptcu_dvnsq /= event_ID_ens.size() ; 
  sum_dpt4_dvnsq /= event_ID_ens.size() ; 
  sum_dpt_dvnsq_dvnsq /= event_ID_ens.size() ; 

  double sigma_pt = sqrt(sum_dptsq);
  double sigma_vnsq = sqrt(sum_dvn2sq);

  double num ; 
  double den ;
   
  num = sum_dpt_dvnsq ; 
  den = pow(sigma_pt,1) * sigma_vnsq ;
  Cov11 = num  ; 
  rho11 = num / den ; 
  
  num = sum_dptsq_dvnsq ; 
  den = pow(sigma_pt,2) * sigma_vnsq ;
  Cov21 = num  ; 
  rho21 = num / den ; 

  num = sum_dptcu_dvnsq - 3 * sum_dptsq * sum_dpt_dvnsq ; 
  den = pow(sigma_pt,3) * sigma_vnsq ;
  Cov31 = num  ; 
  rho31 = num / den ; 

  num = sum_dpt4_dvnsq - 6 * sum_dptsq *  sum_dptsq_dvnsq - 4 * sum_dptcu * sum_dpt_dvnsq ; 
  den = pow(sigma_pt,4) * sigma_vnsq ;
  Cov41 = num  ; 
  rho41 = num / den ; 

  num = sum_dpt_dvnsq_dvnsq - 2 * avg_vnsq * sum_dpt_dvnsq ; 
  den = pow(sigma_pt,1) * pow(sigma_vnsq,2);
  Cov122 = num  ; 
  rho122 = num / den ; 

}



// The observable : Eq(7) of https://arxiv.org/pdf/2506.04029
void observables::output_Bozek_rn_pion_proton(int n){

 // create an ensemble
 std::vector<int> event_ID_ens;
 
 double rn;
 
 double sumC = 0. ; 
 double sumC2 = 0. ; 


 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_Bozek_rn_pion_proton(n,event_ID_ens, rn);

    sumC  += rn ; 
    sumC2 += rn * rn ; 

 }
 
  std::ofstream mFile;
  std::stringstream output_filename;

  // calculate mean and std. dev. of <vn(211)*vnstar(2212)> 
  double Mean_rn = sumC / rmof->get_total_events() ; 
  double Erro_rn = sqrt( sumC2 / rmof->get_total_events() - pow(Mean_rn,2) ) ; 
  // write to file
  output_filename.str("");
  output_filename << "results/Bozek_r" << n << "_pion_proton";
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
  mFile << "#rn_meson_baryon   error " << std::endl ;
  mFile << Mean_rn << "   " << Erro_rn <<  std::endl ; 
  mFile.close();
}


void observables::calculate_Bozek_rn_pion_proton(int n, std::vector<int> event_ID_ens, 
      double& rn){
  // calculate the observable for one ensemble //
  double sum_vn_211_vn_2212_star  = 0. ; 
  double sum_vn_211_vn_211_star   = 0. ; 
  double sum_vn_2212_vn_2212_star = 0. ; 
   
  double temp_vnvnstar;
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 

    temp_vnvnstar = ev->get_integrated_vn(211,n,0) * ev->get_integrated_vn(2212,n,0) 
      + ev->get_integrated_vn(211,n,1) * ev->get_integrated_vn(2212,n,1) ;  
    sum_vn_211_vn_2212_star += temp_vnvnstar ;  

    temp_vnvnstar = ev->get_integrated_vn(211,n,0) * ev->get_integrated_vn(211,n,0) 
      + ev->get_integrated_vn(211,n,1) * ev->get_integrated_vn(211,n,1) ;  
    sum_vn_211_vn_211_star += temp_vnvnstar ;  

    temp_vnvnstar = ev->get_integrated_vn(2212,n,0) * ev->get_integrated_vn(2212,n,0) 
      + ev->get_integrated_vn(2212,n,1) * ev->get_integrated_vn(2212,n,1) ;  
    sum_vn_2212_vn_2212_star += temp_vnvnstar ;  


  }
  
  sum_vn_211_vn_2212_star /= event_ID_ens.size() ; 
  sum_vn_211_vn_211_star /= event_ID_ens.size() ; 
  sum_vn_2212_vn_2212_star /= event_ID_ens.size() ; 

  rn = sum_vn_211_vn_2212_star / sqrt(sum_vn_211_vn_211_star*sum_vn_2212_vn_2212_star) ;                
}





