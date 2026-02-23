#include "event.h"

event::event(int aNptbins):music_pt_bins(aNptbins){

meanpt.resize(NPID, 0.0);

vn.resize(NPID);
for (int i = 0; i < NPID; i++) {
    vn[i].resize(Nharmonics);
    for (int j = 0; j < Nharmonics; j++) {
        vn[i][j].resize(2, 0.0);
    }
}


vnpt.resize(NPID);
for (int pid = 0; pid < NPID; pid++) {
    vnpt[pid].resize(Nharmonics);
    for (int h = 0; h < Nharmonics; h++) {
        vnpt[pid][h].resize(2);
        for (int ri = 0; ri < 2; ri++) {
            vnpt[pid][h][ri].resize(music_pt_bins, 0.0);
        }
    }
}

}
