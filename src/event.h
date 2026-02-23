#pragma once

#include <vector>

class event {

public:

    explicit event(int aNptbins);

    // =========================
    // Getters
    // =========================

    inline double get_integrated_vn(int PID,int h,int ri) const noexcept {
        int pid = get_PID_index(PID);
        return (pid < 0) ? 0.0 : vn[pid][h][ri];
    }

    inline double get_pt_differential_vn(int PID,int h,int ri,int pt) const noexcept {
        int pid = get_PID_index(PID);
        return (pid < 0) ? 0.0 : vnpt[pid][h][ri][pt];
    }

    inline double get_mean_pt(int PID) const noexcept {
        int pid = get_PID_index(PID);
        return (pid < 0) ? 0.0 : meanpt[pid];
    }


    // =========================
    // Setters
    // =========================

    inline void set_integrated_vn(int PID,int h,int ri,double val) noexcept {
        int pid = get_PID_index(PID);
        if(pid >= 0) vn[pid][h][ri] = val;
    }

    inline void set_differential_vn(int PID,int h,int ri,int pt,double val) noexcept {
        int pid = get_PID_index(PID);
        if(pid >= 0) vnpt[pid][h][ri][pt] = val;
    }


    inline void set_mean_pt(int PID,double val) noexcept {
        int pid = get_PID_index(PID);
        if(pid >= 0) meanpt[pid] = val;
    }



    // =========================
    // Size info helpers
    // =========================

    inline int get_Nharmonics() const noexcept { return Nharmonics; }
    inline int get_NPID() const noexcept { return NPID; }
    inline int get_Nptbins() const noexcept { return music_pt_bins; }


private:

    // dimensions
    static constexpr int Nharmonics = 5;
    static constexpr int NPID = 7;

    int music_pt_bins;

    // data containers
    std::vector<std::vector<std::vector<double>>> vn;        // [pid][harmonic][re/im]
    std::vector<double> meanpt;                              // [pid]
    std::vector<std::vector<std::vector<std::vector<double>>>> vnpt; // [pid][harmonic][re/im][pt]


    // =========================
    // PID mapping
    // =========================

    inline int get_PID_index(int PID) const noexcept {
        switch (PID) {
            case 0:      return 0;
            case 211:    return 1;
            case -211:   return 2;
            case 321:    return 3;
            case -321:   return 4;
            case 2212:   return 5;
            case -2212:  return 6;
            default:     return -1;
        }
    }
};

