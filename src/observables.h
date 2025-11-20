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
    void output_meanpt_v2sq_correlation();
    void output_meanpt_v3sq_correlation();
    void output_pt_diff_meanpt_v2v2pt_correlation();
    void output_pt_diff_meanpt_v3v3pt_correlation();

  private :
    std::vector<event*> event_arena;
    read_music_output_files* rmof ; 
    random_gen* rand;
    std::vector<int> get_an_event_ensemble();
    void calculate_meanpt_v2sq_correlation(std::vector<int> , double& );
    void calculate_meanpt_v3sq_correlation(std::vector<int> , double& );
        
    void calculate_pt_diff_meanpt_v2v2pt_correlation(std::vector<int> , std::vector<double>& );
    void calculate_pt_diff_meanpt_v3v3pt_correlation(std::vector<int> , std::vector<double>& );
   
    // kinematics cut     
    int pid ; int yflag ; double rapmin ; 
    double rapmax ; double ptmin ; double ptmax ; 

    
};
