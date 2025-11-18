#include "read_music_output_files.h"
#include "event.h"

read_music_output_files::read_music_output_files(std::vector<std::string>  aa_music_output_paths){
  music_output_paths = aa_music_output_paths ;
     
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
      input_filename << "/vnchpT_y_-0.5_0.5.dat";
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
  std::cout << "pt-bin      pt-val" << std::endl ; 
  for(long unsigned int ii = 0 ; ii < ptval.size(); ii++){
    std::cout << "   " << ii << "       " << ptval[ii] << std::endl ;  
  }
  std::cout << "================" << std::endl ;
  
  // set up all events info
  for(int ii=0; ii<total_music_events; ii++){
    event* ev = new event(music_pt_bins);
    event_arena.push_back(ev); 
  }
  
}


void read_music_output_files::read_pt_integrated_stuff_for_hpm(int yflag, double rapmin, double rapmax, double ptmin, double ptmax){
  std::ifstream file;
  int temp_total_music_events = 0 ;
  double Nch, v1cos, v1sin, v2cos, v2sin, v3cos, v3sin, v4cos, v4sin, dummy ; 
  for(long unsigned int output_path_index=0; output_path_index < music_output_paths.size() ; output_path_index++ ){
    for(int ioutputIDX=0; ioutputIDX < 9999 ; ioutputIDX++ ){ // maximum 10,000 outputs/events in a path should be there
      std::stringstream input_filename;
      input_filename.str(std::string());
      input_filename << music_output_paths[output_path_index].c_str() ;
      input_filename << "/outputs_" << std::setfill('0') << std::setw(4) << ioutputIDX; 
      input_filename << "/vnch_pT_";
      input_filename << ptmin << "_" << ptmax ;
      if(yflag==1){
        input_filename << "_y_" ;
      }
      else{
        input_filename << "_eta_" ;
      }
       input_filename << rapmin << "_" << rapmax ; 
       input_filename << ".dat" ;
      file.open(input_filename.str().c_str(), std::ios::in);
      if(!file){
        continue ; 
      }
      else{
        file.getline(buff,500) ; // header
        file.getline(buff,500) ; // line of interest
        iss = new std::istringstream(buff);
        *iss >> Nch >> v1cos >> v1sin >> v2cos >> v2sin >> v3cos >> v3sin >> v4cos >> v4sin >> dummy ;
        event_arena[temp_total_music_events]->set_integrated_vn(0,0,Nch) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(0,1,0) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(1,0,v1cos) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(1,1,v1sin) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(2,0,v2cos) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(2,1,v2sin) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(3,0,v3cos) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(3,1,v3sin) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(4,0,v4cos) ; 
        event_arena[temp_total_music_events]->set_integrated_vn(4,1,v4sin) ; 
        delete iss ; 
        file.close();
        temp_total_music_events ++ ; 
      }
    } // iouputIdx
  } // loop over paths
  
  if(temp_total_music_events == 0){
    std::cout << "reading pt integrated quantities failure !!! Not got a single file of the given kinematic cut" << std::endl ; 
    exit(-1);
  }
  
    if(temp_total_music_events > total_music_events || temp_total_music_events < total_music_events){
     std::cout << "reading pt integrated quantities failure !!! unexpectedly reading more/less events" << std::endl ; 
     exit(-1);
   }
  
}



void read_music_output_files::read_meanpt_for_hpm(int yflag, double rapmin, double rapmax, double ptmin, double ptmax){
  std::ifstream file;
  int temp_total_music_events = 0 ;
  double mpt ; 
  for(long unsigned int output_path_index=0; output_path_index < music_output_paths.size() ; output_path_index++ ){
    for(int ioutputIDX=0; ioutputIDX < 9999 ; ioutputIDX++ ){ // maximum 10,000 outputs/events in a path should be there
      std::stringstream input_filename;
      input_filename.str(std::string());
      input_filename << music_output_paths[output_path_index].c_str() ;
      input_filename << "/outputs_" << std::setfill('0') << std::setw(4) << ioutputIDX; 
      input_filename << "/mean_pt_ch_pt_";
      input_filename << ptmin << "_" << ptmax ;
      if(yflag==1){
        input_filename << "_y_" ;
      }
      else{
        input_filename << "_eta_" ;
      }
       input_filename << rapmin << "_" << rapmax ; 
       input_filename << ".dat" ;
      file.open(input_filename.str().c_str(), std::ios::in);
      if(!file){
        continue ; 
      }
      else{
        file.getline(buff,500) ; // header
        file.getline(buff,500) ; // line of interest
        iss = new std::istringstream(buff);
        *iss >> mpt ;
        event_arena[temp_total_music_events]->set_meanpt(mpt) ; 
        delete iss ; 
        file.close();
        temp_total_music_events ++ ; 
      }
    } // iouputIdx
  } // loop over paths
  
    if(temp_total_music_events == 0){
     std::cout << "reading mean-pt failure !!! Not got a single file of the given kinematic cut" << std::endl ; 
     exit(-1);
    }
    if(temp_total_music_events > total_music_events || temp_total_music_events < total_music_events){
     std::cout << "reading mean-pt failure !!! unexpectedly reading more/less events" << std::endl ; 
     exit(-1);
   }
  
}



void read_music_output_files::read_pt_differential_stuff_for_hpm(int yflag, double rapmin, double rapmax){
  std::ifstream file;
  int temp_total_music_events = 0 ;
  double ptv, dnptdptdy, v1cos, v1sin, v2cos, v2sin, v3cos, v3sin, v4cos, v4sin, dummy ; 
  for(long unsigned int output_path_index=0; output_path_index < music_output_paths.size() ; output_path_index++ ){
    for(int ioutputIDX=0; ioutputIDX < 9999 ; ioutputIDX++ ){ // maximum 10,000 outputs/events in a path should be there
      std::stringstream input_filename;
      input_filename.str(std::string());
      input_filename << music_output_paths[output_path_index].c_str() ;
      input_filename << "/outputs_" << std::setfill('0') << std::setw(4) << ioutputIDX; 
      input_filename << "/vnchpT";
      if(yflag==1){
        input_filename << "_y_" ;
      }
      else{
        input_filename << "_eta_" ;
      }
       input_filename << rapmin << "_" << rapmax ; 
       input_filename << ".dat" ;
      file.open(input_filename.str().c_str(), std::ios::in);
      if(!file){
        continue ; 
      }
      else{
        file.getline(buff,500) ; // header
        int ii=0;
        while (ii < 100 &&  file.getline(buff,500)){           
          iss = new std::istringstream(buff);
          *iss >> ptv >> dnptdptdy >> v1cos >> v1sin >> v2cos >> v2sin >> v3cos >> v3sin >> v4cos >> v4sin >> dummy ;
          if(temp_total_music_events==1  && ii==5){ 
            if( fabs(ptv-ptval[ii]) > 0.0001){
               std::cout << "pt bin error ...   ii = "  
               << ii << "  pt = " << ptval[ii] 
               << "  readptv = " << ptv << "   " 
               << input_filename.str() 
               << std::endl; 
               exit(-1);
             }   
          }
          event_arena[temp_total_music_events]->set_differential_vn(0,0,ii,dnptdptdy) ; 
          event_arena[temp_total_music_events]->set_differential_vn(0,1,ii,0) ; 
          event_arena[temp_total_music_events]->set_differential_vn(1,0,ii,v1cos) ; 
          event_arena[temp_total_music_events]->set_differential_vn(1,1,ii,v1sin) ; 
          event_arena[temp_total_music_events]->set_differential_vn(2,0,ii,v2cos) ; 
          event_arena[temp_total_music_events]->set_differential_vn(2,1,ii,v2sin) ; 
          event_arena[temp_total_music_events]->set_differential_vn(3,0,ii,v3cos) ; 
          event_arena[temp_total_music_events]->set_differential_vn(3,1,ii,v3sin) ; 
          event_arena[temp_total_music_events]->set_differential_vn(4,0,ii,v4cos) ; 
          event_arena[temp_total_music_events]->set_differential_vn(4,1,ii,v4sin) ; 
          delete iss ;
          ii++ ; 
        } 
        file.close();
        temp_total_music_events ++ ; 
      }
    } // iouputIdx
  } // loop over paths
  
  if(temp_total_music_events == 0){
    std::cout << "reading pt differential failure !!! Not got a single file of the given kinematic cut" << std::endl ; 
    exit(-1);
  }
  if(temp_total_music_events > total_music_events || temp_total_music_events < total_music_events ){
     std::cout << "reading pt differential failure !!! unexpectedly reading more/less events" << std::endl ; 
     exit(-1);
   }
}




