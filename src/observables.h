#pragma once
#include "event.h"
#include "read_music_output_files.h"

class observables{

  public :

    observables(read_music_output_files* );
    void calculate_meanpt_v2sq_correlation();
 
  private :
        std::vector<event*> event_arena;
        read_music_output_files* rmof ; 
    
};
