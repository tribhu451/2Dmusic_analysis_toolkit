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

    observables(read_music_output_files*, int yflag, double rapmin, double rapmax, double ptmin, double ptmax);
    void output_meanpt_vnsq_correlation_charged_hadrons(int n);
    void output_pt_diff_meanpt_vnvnpt_correlation_charged_hadrons(int n);
    void output_pt_diff_meanpt_vnptvnpt_correlation_charged_hadrons(int n);
    void output_pt_diff_multiparticle_vn_charged_hadrons(int n);
    void output_pt_diff_multiparticle_vn_method2_charged_hadrons(int n);
    void output_average_multiplicity_charged_hadrons();
    void output_meanpt_vnsq_higher_moments_charged_hadrons(int n);
    void output_meanpt_vnsq_higher_moments_mult_fluc_corrected_charged_hadrons(int n);
    void output_Bozek_rn_pion_proton(int n);


  private :
    std::vector<event*> event_arena;
    read_music_output_files* rmof ; 
    random_gen* rand;
    std::vector<int> get_an_event_ensemble();
    void calculate_meanpt_vnsq_correlation_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      double&  Mpt, double&  M_ptsq,  double&  M_vnvnstar,  double&  M_vnvnstarsq, 
      double& Cov_Mpt_vnvnstar,   double& meanpt_vnsq_corr_of_one_ens);          
    void calculate_pt_diff_meanpt_vnvnpt_correlation_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      std::vector<double>& M_vnvnptstar, std::vector<double>& M_vnvnptstar_sq, 
      std::vector<double>& cov, std::vector<double>& obs);
    void calculate_pt_diff_meanpt_vnptvnpt_correlation_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      std::vector<double>& M_vnptvnptstar, std::vector<double>& M_vnptvnptstar_sq, 
      std::vector<double>& cov, std::vector<double>& obs);
    void calculate_pt_diff_multiparticle_vn_charged_hadrons(int n, std::vector<int> event_ID_ens, double&, double&, 
      std::vector<double>& vn_2, std::vector<double>& vn_4);
    void calculate_pt_diff_multiparticle_vn_method2_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      double& vn_sq, double& vn_fr, std::vector<double>& vn_2_num, std::vector<double>& vn_4_num );
    void calculate_meanpt_vnsq_higher_moments_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      double& Cov11, double& rho11, double& Cov21, double& rho21, 
      double& Cov31, double& rho31, double& Cov41, double& rho41, 
      double& Cov122, double& rho122 );
    void calculate_meanpt_vnsq_higher_moments_mult_fluc_corrected_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      double& Cov11, double& rho11, double& Cov21, double& rho21, 
      double& Cov31, double& rho31, double& Cov41, double& rho41, 
      double& Cov122, double& rho122 );
    void calculate_Bozek_rn_pion_proton(int n, std::vector<int> event_ID_ens, 
      double& rn);
   
    // kinematics cut     
    int yflag ; double rapmin ; 
    double rapmax ; double ptmin ; double ptmax ; 

    
};
