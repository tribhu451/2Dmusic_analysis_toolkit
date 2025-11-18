#pragma once
#include "event.h"
#include "read_music_output_files.h"
#include "random.h"
#include <vector>

class observables{

  public :

    observables(read_music_output_files* );
    void output_meanpt_v2sq_correlation();
    void output_meanpt_v3sq_correlation();
 
  private :
        std::vector<event*> event_arena;
        read_music_output_files* rmof ; 
        random_gen* rand;
        std::vector<int> get_an_event_ensemble();
        void calculate_meanpt_v2sq_correlation(std::vector<int> , double& );
        void calculate_meanpt_v3sq_correlation(std::vector<int> , double& );
    
};
