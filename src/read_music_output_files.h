#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "event.h"

class read_music_output_files{

  public :
    read_music_output_files(std::vector<std::string>  aa, int pid, int yflag, double rapmin, double rapmax, double ptmin, double ptmax);
    inline int get_music_pit_bins(){return music_pt_bins;};
    inline int get_total_music_events(){return total_music_events;};
    std::vector<double> get_music_ptvals(){return ptval;};
    void read_pt_integrated_stuff();
    void read_meanpt();
    void read_pt_differential_stuff();
    std::vector<event*> get_event_arena(){return event_arena;};
    event* get_event(int ii){return event_arena[ii];}
    double get_pt_val_of_bin(int ii){return ptval[ii];}
    
  private :
    std::istringstream* iss;
    char buff[400];
    std::vector<std::string> music_output_paths ; 
    int music_pt_bins ; 
    int total_music_events ; 
    std::vector<double> ptval ; 
    std::vector<event*> event_arena;
    int pid ; int yflag ; double rapmin ; 
    double rapmax ; double ptmin ; double ptmax ; 
};
