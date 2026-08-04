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











void observables::output_vo_ratio_proton_hpm(){
   std::vector<int> event_ID_ens;
   std::ofstream mFile;
   std::stringstream output_filename;
   double voratio ;
   double sumx = 0. ; 
   double sumx2 = 0. ; 
   for(int ii=0; ii<rmof->get_total_events(); ii++){
     event_ID_ens = get_an_event_ensemble();
     // calculate correlation of a given ensemble
     calculate_vo_ratio_proton_hpm(event_ID_ens, voratio);
     sumx  += voratio ; 
     sumx2 += pow(voratio,2); 
   }

   sumx  /= rmof->get_total_events();
   sumx2 /= rmof->get_total_events();
 
   output_filename.str("");
   output_filename << "results/";
   output_filename << "vo_ratio_proton_";
   output_filename << "charged_hadrons";
   output_filename << "_pt_";
   output_filename << ptmin << "_" << ptmax ;
   if(yflag==1){
    output_filename << "_y_" ; }
   else{
    output_filename << "_eta_" ;}
   output_filename << rapmin << "_" << rapmax ;
   output_filename << ".dat";
   mFile.open(output_filename.str().c_str(), std::ios::out );

   mFile << "#v0_ratio(p/ch)  error" << std::endl ;
   mFile << sumx << "  " << sqrt( sumx2 - sumx * sumx ) << std::endl ; 
   mFile.close(); 
}


void observables::calculate_vo_ratio_proton_hpm( std::vector<int> event_ID_ens, 
     double&  voratio){
  // calculate the observable for one ensemble //
  double mpt_ch = 0. ; 
  double mpt_prot = 0. ; 

  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    mpt_ch += ev->get_mean_pt(0);
    mpt_prot += ev->get_mean_pt(2212);
  }
  
  mpt_ch     /= event_ID_ens.size() ; 
  mpt_prot   /= event_ID_ens.size() ; 

  double deltaptsq_ch = 0 ;
  double deltaptsq_prot = 0 ; 
  double spt ; 
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    spt = ev->get_mean_pt(0);
    deltaptsq_ch += pow(spt - mpt_ch,2) ; 
    spt = ev->get_mean_pt(2212);
    deltaptsq_prot += pow(spt - mpt_prot,2) ; 
  }

  deltaptsq_ch     /= event_ID_ens.size() ; 
  deltaptsq_prot   /= event_ID_ens.size() ; 

  voratio = ( sqrt(deltaptsq_prot) / mpt_prot ) 
   /  ( sqrt(deltaptsq_ch) / mpt_ch ) ; 

}



void observables::output_ebe_meanpt_correlation_hpm_proton(){
   std::vector<int> event_ID_ens;
   std::ofstream mFile;
   std::stringstream output_filename;
   double pearson;
   double sumx = 0. ; 
   double sumx2 = 0. ; 
   for(int ii=0; ii<rmof->get_total_events(); ii++){
     event_ID_ens = get_an_event_ensemble();
     // calculate correlation of a given ensemble
     calculate_ebe_meanpt_correlation_hpm_proton(event_ID_ens,pearson);
     sumx  += pearson ; 
     sumx2 += pow(pearson,2); 
   }

    sumx  /= rmof->get_total_events();
    sumx2 /= rmof->get_total_events();

    output_filename.str("");
    output_filename << "results/";
    output_filename << "ebe_meanpt_pearson_correlation_hpm_proton";
    output_filename << "_pt_";
    output_filename << ptmin << "_" << ptmax ;
    if(yflag==1){
     output_filename << "_y_" ; }
    else{
     output_filename << "_eta_" ;}
    output_filename << rapmin << "_" << rapmax ;
    output_filename << ".dat";
    mFile.open(output_filename.str().c_str(), std::ios::out );

    mFile << "#rho([pt](pi+),[pt](p))    error  " << std::endl ;
    mFile << sumx << "  " << sqrt( sumx2 - sumx * sumx ) << std::endl ; 
    mFile.close();
  
}

void observables::calculate_ebe_meanpt_correlation_hpm_proton(std::vector<int> event_ID_ens, 
     double&  pearson){
  // calculate the observable for one ensemble //
  double sumxy = 0. ; 
  double sumx  = 0. ; 
  double sumx2 = 0. ; 
  double sumy  = 0. ; 
  double sumy2 = 0. ;
      
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    double ptHpm = ev->get_mean_pt(0);
    double ptProt = ev->get_mean_pt(2212);
    sumx  += ptHpm ; 
    sumx2 += ptHpm * ptHpm ; 
    sumy  += ptProt ; 
    sumy2 += ptProt * ptProt ; 
    sumxy += ptHpm * ptProt ; 
  }
  
  sumx /= event_ID_ens.size() ; 
  sumx2 /= event_ID_ens.size() ;
  sumy /= event_ID_ens.size() ; 
  sumy2 /= event_ID_ens.size() ;
  sumxy /= event_ID_ens.size() ; 

  pearson =  (sumxy - sumx * sumy) / ( sqrt(sumx2 - sumx * sumx) * sqrt(sumy2 - sumy * sumy) ) ; 
}





// The observable : Eq(7) of https://arxiv.org/pdf/2506.04029
void observables::output_Bozek_rn(int n, int PID1, int PID2){

 // create an ensemble
 std::vector<int> event_ID_ens;
 
 double rn;
 
 double sumC = 0. ; 
 double sumC2 = 0. ; 


 for(int ii=0; ii<rmof->get_total_events(); ii++){
    event_ID_ens = get_an_event_ensemble();
    // calculate correlation of a given ensemble
    calculate_Bozek_rn(n,PID1, PID2,event_ID_ens, rn);

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
  output_filename << "results/Bozek_r" << n << "_" ;
  if(PID1==0)
     output_filename << "hpm_" ;
  else
    output_filename << PID1 << "_" ; 
  if(PID2==0)
     output_filename << "hpm" ;
  else
    output_filename << PID2 ; 
    
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
  mFile << "#rn_hpm_proton   error " << std::endl ;
  mFile << Mean_rn << "   " << Erro_rn <<  std::endl ; 
  mFile.close();
}




void observables::calculate_Bozek_rn(int n, int PID1, int PID2, std::vector<int> event_ID_ens, 
      double& rn){
  // calculate the observable for one ensemble //
  double sum_vn_PID1_vn_PID2_star  = 0. ; 
  double sum_vn_PID1_vn_PID1_star   = 0. ; 
  double sum_vn_PID2_vn_PID2_star = 0. ; 
   
  double temp_vnvnstar;
   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 

    temp_vnvnstar = ev->get_integrated_vn(PID1,n,0) * ev->get_integrated_vn(PID2,n,0) 
      + ev->get_integrated_vn(PID1,n,1) * ev->get_integrated_vn(PID2,n,1) ;  
    sum_vn_PID1_vn_PID2_star += temp_vnvnstar ;  

    temp_vnvnstar = ev->get_integrated_vn(PID1,n,0) * ev->get_integrated_vn(PID1,n,0) 
      + ev->get_integrated_vn(PID1,n,1) * ev->get_integrated_vn(PID1,n,1) ;  
    sum_vn_PID1_vn_PID1_star += temp_vnvnstar ;  

    temp_vnvnstar = ev->get_integrated_vn(PID2,n,0) * ev->get_integrated_vn(PID2,n,0) 
      + ev->get_integrated_vn(PID2,n,1) * ev->get_integrated_vn(PID2,n,1) ;  
    sum_vn_PID2_vn_PID2_star += temp_vnvnstar ;  


  }
  
  sum_vn_PID1_vn_PID2_star /= event_ID_ens.size() ; 
  sum_vn_PID1_vn_PID1_star /= event_ID_ens.size() ; 
  sum_vn_PID2_vn_PID2_star /= event_ID_ens.size() ; 

  rn = sum_vn_PID1_vn_PID2_star / sqrt(sum_vn_PID1_vn_PID1_star*sum_vn_PID2_vn_PID2_star) ;                
}




void observables::output_vn_meanpt_correlation_ratio_proton_charged_hadron(int n){
   std::vector<int> event_ID_ens;
   std::ofstream mFile;
   std::stringstream output_filename;
   double rhoCh, rhoP, ratio ;
   double sumu = 0. ; 
   double sumu2 = 0. ; 
   double sumv = 0. ; 
   double sumv2 = 0. ; 
   double sumx = 0. ; 
   double sumx2 = 0. ; 
   for(int ii=0; ii<rmof->get_total_events(); ii++){
     event_ID_ens = get_an_event_ensemble();
     // calculate correlation of a given ensemble
     calculate_vn_meanpt_correlation_ratio_proton_charged_hadron(event_ID_ens, n, rhoCh, rhoP, ratio);
     sumu  += rhoCh ; 
     sumu2 += pow(rhoCh,2); 
     sumv  += rhoP ; 
     sumv2 += pow(rhoP,2); 
     sumx  += ratio ; 
     sumx2 += pow(ratio,2); 
   }

   sumu  /= rmof->get_total_events();
   sumu2 /= rmof->get_total_events();
   sumv  /= rmof->get_total_events();
   sumv2 /= rmof->get_total_events();
   sumx  /= rmof->get_total_events();
   sumx2 /= rmof->get_total_events();
  
   output_filename.str("");
   output_filename << "results/";
   output_filename << "v" << n <<"_meanpt_correlation_ratio_proton_";
   output_filename << "charged_hadrons";
   output_filename << "_pt_";
   output_filename << ptmin << "_" << ptmax ;
   if(yflag==1){
    output_filename << "_y_" ; }
   else{
    output_filename << "_eta_" ;}
   output_filename << rapmin << "_" << rapmax ;
   output_filename << ".dat";
   mFile.open(output_filename.str().c_str(), std::ios::out );

   mFile << "#rho_ch    error    rho_proton    error    vn_meanpt_correlation_ratio(p/ch)  error" << std::endl ;
   mFile << sumu << "  " << sqrt( sumu2 - sumu * sumu ) << "  " 
     << sumv << "  " << sqrt( sumv2 - sumv * sumv )  << "  "
     << sumx << "  " << sqrt( sumx2 - sumx * sumx )  << std::endl ; 
   mFile.close(); 
}


void observables::calculate_vn_meanpt_correlation_ratio_proton_charged_hadron
(std::vector<int> event_ID_ens, int n,  double&  rhoch,  double&  rhop, double&  RatioRho){
  calculate_vn_meanpt_correlation(0,n,event_ID_ens,rhoch);
  calculate_vn_meanpt_correlation(2212,n,event_ID_ens,rhop);
  RatioRho = rhop/rhoch ;
}


void observables::calculate_vn_meanpt_correlation( int PID, int n, std::vector<int> event_ID_ens, 
     double&  rho){

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
    nch = ev->get_integrated_vn(0,0,0) ;
    spt = ev->get_mean_pt(PID);
    vnvnstar = pow(ev->get_integrated_vn(PID,n,0) ,2) + pow(ev->get_integrated_vn(PID,n,1),2) ;  
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
  double avg_vnsq = sumvnvnstar ;                                // <vn^2> charged hadron

  // second : cov(delta pt, Nch) and cov(delta vn^2, Nch) //
  double sum_deltapt_nch = 0 ; 
  double sum_deltavnsq_nch = 0 ; 
  double sum_deltapt = 0 ; 
  double sum_deltavnsq = 0 ;
  double sum_nch = 0 ;  
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    nch =  ev->get_integrated_vn(0,0,0) ;
    spt = ev->get_mean_pt(PID);
    vnvnstar = pow(ev->get_integrated_vn(PID,n,0) ,2) + pow(ev->get_integrated_vn(PID,n,1),2) ;  
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



  double num ; 
  double den ;
  double sum_dptsq; 
  double sum_dpt_dvnsq; 
  double deltapt_tilde, delta_vnsq_tilde;
  
  // multiplicity fluctuation corrected rho([pt],v2) calculation for charged hadron //
  sum_dptsq=0.; 
  sum_dpt_dvnsq=0.;   
  for(long unsigned int ii=0; ii<event_ID_ens.size(); ii++){
    int eventID = event_ID_ens[ii] ; 
    event* ev = rmof->get_event(eventID) ; 
    spt = ev->get_mean_pt(PID);
    nch = ev->get_integrated_vn(0,0,0) ;
    vnvnstar = pow(ev->get_integrated_vn(PID,n,0),2) + pow(ev->get_integrated_vn(PID,n,1),2) ;  
    deltapt_tilde =  (spt-avg_mpt) - cov_deltapt_nch / variance_nch * ( nch - avg_nch ) ; 
    delta_vnsq_tilde =  (vnvnstar-avg_vnsq) - cov_deltavnsq_nch / variance_nch * ( nch - avg_nch ) ; 
    sum_dptsq += pow(deltapt_tilde,2.);
    sum_dpt_dvnsq +=  (deltapt_tilde)*(delta_vnsq_tilde);
  }
  
  sum_dptsq /= event_ID_ens.size() ; 
  sum_dpt_dvnsq /= event_ID_ens.size() ; 
  
  double sigma_pt = sqrt(sum_dptsq);

  num = sum_dpt_dvnsq ; 
  den = sigma_pt * avg_vnsq ;
  rho = num / den ;   // rho([pt],v2) 
}
























