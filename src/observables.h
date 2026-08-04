#pragma once
#include "event.h"
#include "read_music_output_files.h"
#include "random.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

class observables{

  public :

    observables(read_music_output_files* , int pid, int yflag, double rapmin, double rapmax, double ptmin, double ptmax);
    void output_meanpt_vnsq_correlation(int n);
    void output_pt_diff_meanpt_vnvnpt_correlation(int n);
    void output_pt_diff_meanpt_vnptvnpt_correlation(int n);
    void output_pt_diff_multiparticle_vn(int n);
    void output_pt_diff_multiparticle_vn_method2(int n);
    void output_meanpt_vnsq_higher_moments_mult_fluc_corrected(int n);
    void output_v0_vn_in_both_w_and_wo_mult_fluc_correction(int n);
    void output_relation(int n);
    void output_moments_vn(int n);   
     
  private :
    std::vector<event*> event_arena;
    read_music_output_files* rmof ; 
    random_gen* rand;
    std::vector<int> get_an_event_ensemble();
    void calculate_meanpt_vnsq_correlation(int n, std::vector<int> event_ID_ens, 
      double&  Mpt, double&  M_ptsq,  double&  M_vnvnstar,  double&  M_vnvnstarsq, 
      double& Cov_Mpt_vnvnstar,   double& meanpt_vnsq_corr_of_one_ens);          
    void calculate_pt_diff_meanpt_vnvnpt_correlation(int n, std::vector<int> event_ID_ens, 
      std::vector<double>& M_vnvnptstar, std::vector<double>& M_vnvnptstar_sq, 
      std::vector<double>& cov, std::vector<double>& obs);
    void calculate_pt_diff_meanpt_vnptvnpt_correlation(int n, std::vector<int> event_ID_ens, 
      std::vector<double>& M_vnptvnptstar, std::vector<double>& M_vnptvnptstar_sq, 
      std::vector<double>& cov, std::vector<double>& obs);
    void calculate_pt_diff_multiparticle_vn(int n, std::vector<int> event_ID_ens, double&, double&, 
      std::vector<double>& vn_2, std::vector<double>& vn_4);
    void calculate_pt_diff_multiparticle_vn_method2(int n, std::vector<int> event_ID_ens, 
      double& vn_sq, double& vn_fr, std::vector<double>& vn_2_num, std::vector<double>& vn_4_num );
    void calculate_meanpt_vnsq_higher_moments_mult_fluc_corrected(int n, std::vector<int> event_ID_ens, 
      double& Cov11, double& rho11, double& Cov21, double& rho21, 
      double& Cov31, double& rho31, double& Cov41, double& rho41, 
      double& Cov12, double& rho12, double& Cov13, double& rho13 );
    void calculate_v0_vn_in_both_w_and_wo_mult_fluc_correction(int n, std::vector<int> event_ID_ens, 
    double& Mpt_wo, double& Sigmapt_wo, double& v0_wo, double& Mpt_w, double& Sigmapt_w, double& v0_w,
    double& Mvnsq_wo, double& Sigmavnsq_wo, double& ttvn_wo, double& Mvnsq_w, double& Sigmavnsq_w, double& ttvn_w );
    void calculate_relation(int n, std::vector<int> event_ID_ens, double& , double&  );
    void calculate_moments_vn(int n, std::vector<int> event_ID_ens, double& , double&, double&  );

    // kinematics cut     
    int pid ; int yflag ; double rapmin ; 
    double rapmax ; double ptmin ; double ptmax ; 

    
};
