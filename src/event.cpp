#include "event.h"

event::event(int aNptbins):meanpt(0.),music_pt_bins(aNptbins){
  for(int ii=0; ii<Nharmonics; ii++){
    for(int jj=0; jj<2; jj++){
      vn[ii][jj] = 0. ; 
    }
  }
  
  vnpt.resize(Nharmonics); // first dimension
  for (int ii = 0; ii < Nharmonics; ii++) {
    vnpt[ii].resize(2); // second dimension 
  }
  
  for (int ii = 0; ii < music_pt_bins; ii++) { //  third dimension
    for (int jj = 0; jj < Nharmonics; jj++) {
      for (int kk = 0; kk < 2; kk++) {
        vnpt[jj][kk].push_back(0.0);
      }
    }
  }
  
  
}
