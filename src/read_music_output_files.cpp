#include "read_music_output_files.h"
#include "event.h"

read_music_output_files::read_music_output_files(std::vector<std::string>  aa_music_output_paths, int mNpt): music_pt_bins(mNpt){
  music_output_paths = aa_music_output_paths ;
  
  event* ev = new event(music_pt_bins); 
   
  // now set how many events do you have
  // and also set Number of ptbins ans pt-values ?
  std::ifstream file;
  int temp_total_music_events = 0 ;
  int temp_pt_bins=0 ; 
  for(long unsigned int output_path_index=0; output_path_index < music_output_paths.size() ; output_path_index++ ){
    for(int ioutputIDX=0; ioutputIDX < 9999 ; ioutputIDX++ ){ // maximum 10,000 outputs/events in a path should be there
      std::stringstream input_filename;
      input_filename.str(std::string());
      input_filename << music_output_paths[output_path_index].c_str() ;
      input_filename << "/outputs_" << std::setfill('0') << std::setw(4) << ioutputIDX; 
      input_filename << "/Fvnpt_y-211.dat";
      file.open(input_filename.str().c_str(), std::ios::in);
      if(!file){
        continue ; 
      }
      else{
        temp_total_music_events ++ ; 
        if(temp_total_music_events==1){
          file.getline(buff,450) ; // header
          int ii=0; 
          double tempptval, dummy ; 
          while (ii < 100 &&  file.getline(buff,500)){           
            iss = new std::istringstream(buff);
            *iss >> tempptval >> dummy ;
            //std::cout << "pt = " << tempptval << std::endl ; 
            ptval.push_back(tempptval) ; 
            delete iss;
            temp_pt_bins++ ; 
            ii++;
          }
        }
        file.close();
      }
    }
  }
  total_music_events = temp_total_music_events ; 
  music_pt_bins = temp_pt_bins ; 
  std::cout<< "total MUSIC events = " << total_music_events  << std::endl ; 
  std::cout<< "total MUSIC pt bins = " << ptval.size()  << std::endl ; 
  for(int ii = 0 ; ii < ptval.size(); ii++){
    std::cout <<  "ptbin = " << ii << "  ptval = " << ptval[ii] << std::endl ;  
  }
}













