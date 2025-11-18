#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <iomanip>
#include <sstream>

class read_music_output_files{

  public :
    read_music_output_files(std::vector<std::string>  aa, int);
    inline int get_music_pit_bins(){return music_pt_bins;};
    inline int get_total_music_events(){return total_music_events;};
    std::vector<double> get_music_ptvals(){return ptval;};
 
  private :
    std::istringstream* iss;
    char buff[400];
    std::vector<std::string> music_output_paths ; 
    int music_pt_bins ; 
    int total_music_events ; 
    std::vector<double> ptval ; 
};
