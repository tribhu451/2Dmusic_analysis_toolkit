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

    void output_dndy_or_dndeta();
    void output_pt_diff_multiparticle_vn_charged_hadrons(int n);
    void output_pt_diff_multiparticle_vn_method2_charged_hadrons(int n);

    void output_vo_ratio_proton_hpm();
    void output_ebe_meanpt_correlation_hpm_proton();
    void output_Bozek_rn(int n,int PID1, int PID2);
    void output_vn_meanpt_correlation_ratio_proton_charged_hadron(int n);


  private :
    std::vector<event*> event_arena;
    read_music_output_files* rmof ; 
    random_gen* rand;
    std::vector<int> get_an_event_ensemble();
    void calculate_pt_diff_multiparticle_vn_charged_hadrons(int n, std::vector<int> event_ID_ens, double&, double&, 
      std::vector<double>& vn_2, std::vector<double>& vn_4);
    void calculate_pt_diff_multiparticle_vn_method2_charged_hadrons(int n, std::vector<int> event_ID_ens, 
      double& vn_sq, double& vn_fr, std::vector<double>& vn_2_num, std::vector<double>& vn_4_num );
    void calculate_Bozek_rn(int n, int PID1, int PID2, std::vector<int> event_ID_ens, 
      double& rn);
    void calculate_vo_ratio_proton_hpm( std::vector<int> event_ID_ens, double&  voratio);
    void calculate_ebe_meanpt_correlation_hpm_proton(std::vector<int> event_ID_ens, double& person); 
    
    
    void calculate_vn_meanpt_correlation( int PID, int n, std::vector<int> event_ID_ens, double&  rho);
    void calculate_vn_meanpt_correlation_ratio_proton_charged_hadron
      (std::vector<int> event_ID_ens, int n,  double&  rhoch,  double&  rhop, double&  RatioRho);   
    // kinematics cut     
    int yflag ; double rapmin ; 
    double rapmax ; double ptmin ; double ptmax ; 

    
};
