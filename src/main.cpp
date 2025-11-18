#include<iostream>
#include<fstream>
#include <vector>
#include "event.h"
#include "read_music_output_files.h"
#include "observables.h"


int main(int argc, char **argv){

 int MUSICptbins = 30 ; 
 int yflag = 1 ; 
 double rapmin = -0.5 ; 
 double rapmax =  0.5 ; 
 double ptmin  =  0.2 ; 
 double ptmax  =  3.0 ; 

 std::cout << "========================" << std::endl ; 
 std::cout << " MUSIC analysis-toolkit" << std::endl ; 
 std::cout << "========================" << std::endl ; 
 std::cout << "\n" << std::endl ; 

 if(argc < 2){
    std::cout << "atleast one path required ..." << std::endl ;
    exit(-1); 
  }

 if(argc > 6){
    std::cout << "too many paths provided (5 allowed) ..." << std::endl ;
    exit(-1); 
  }


 std::vector<std::string> music_output_paths ; 
 
 for(int ii=2; ii<7; ii++){
   if(ii==argc){
     std::cout <<  (ii-1) << " paths provided" << std::endl ; 
     for(int jj=1; jj<argc; jj++){
       music_output_paths.push_back(argv[jj]) ;
     } 
     break ; 
   }
 }
 
 
 std::cout << "path names : " ; 
 for(long unsigned int ii=0; ii<music_output_paths.size(); ii++){
   std::cout << music_output_paths[ii] << "/,  " ; 
 }
 std::cout << std::endl ;
 
 
 read_music_output_files* rmof = new read_music_output_files(music_output_paths,MUSICptbins); 
 rmof->read_pt_differential_stuff_for_hpm(yflag, rapmin, rapmax);
 rmof->read_pt_integrated_stuff_for_hpm( yflag,  rapmin,  rapmax,  ptmin,  ptmax );
 rmof->read_meanpt_for_hpm( yflag,  rapmin,  rapmax,  ptmin,  ptmax );
 
 observables* obj = new observables(rmof);
 obj->output_meanpt_v2sq_correlation();
 obj->output_meanpt_v3sq_correlation();
 
 return 0 ;  
}




