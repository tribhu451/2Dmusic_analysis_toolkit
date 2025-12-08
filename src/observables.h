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

    // kinematics cut     
    int pid ; int yflag ; double rapmin ; 
    double rapmax ; double ptmin ; double ptmax ; 

    
};
