#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>

class event{

  public :

    event(int );
    inline double get_integrated_vn(int iharmonics, int real_img ){return vn[iharmonics][real_img];}
    inline double get_pt_differential_vn(int iharmonics, int real_img, int iptbin ){return vnpt[iharmonics][real_img][iptbin];}
    inline double get_mean_pt(){return meanpt;};
 
    inline void set_integrated_vn(int iharmonics, int real_img, double val ){ vn[iharmonics][real_img] = val;}
    inline void set_differential_vn(int iharmonics, int real_img, int ptbinIDX, double val ){ vnpt[iharmonics][real_img][ptbinIDX] = val;}
    inline void set_meanpt(double val){meanpt=val;}
    
  private :
    int Nharmonics = 5 ; 
    
    // integrated quantities
    double vn[5][2] ; // first index is harmonic number where 0 is for yield, second index is real or imaginary part
    double meanpt ;
  
    // pt differential quantities
    int music_pt_bins;
    std::vector<std::vector<std::vector<double>>> vnpt;
    
};
