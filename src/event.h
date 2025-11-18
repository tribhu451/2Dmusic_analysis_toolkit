#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>

class event{

  public :

    event(int );
 
  private :
    int Nharmonics = 5 ; 
    
    // integrated quantities
    double vn[5][2] ; // first index is harmonic number where 0 is for yield, second index is real or imaginary part
    double meanpt ;
  
    // pt differential quantities
    int music_pt_bins;
    std::vector<std::vector<std::vector<double>>> vnpt;
    
};
