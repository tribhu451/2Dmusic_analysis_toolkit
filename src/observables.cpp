#include "observables.h"

observables::observables(read_music_output_files* armof) : rmof(armof){
  event_arena = rmof->get_event_arena();
  
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


void observables::calculate_meanpt_v2sq_correlation(){

 // create an ensemble
 
 
 // calculate correlation of a given ensemble
 
 
 // calculate mean and std. dev.
 
 
 // print the output

}

