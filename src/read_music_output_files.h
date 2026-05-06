#pragma once

#include <vector>
#include <string>
#include "event.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>


class read_music_output_files {

public:

    read_music_output_files(
        const std::vector<std::string>& paths,
        int max_nevents,
        int ayflag,
        double arapmin,double arapmax,
        double aptmin,double aptmax
    );

    ~read_music_output_files();

    void compute_differential_vn_charged_hadron();

    // =============================
    // calculators (return values)
    // =============================
    double calc_meanpt(event* ev,int PID,
                       double ptmin,double ptmax);

    double calc_integrated_vn(event* ev,int PID,
                              int harmonic,int ri,
                              double ptmin,double ptmax);

    // =============================
    // setters for species
    // ==============================
    void compute_meanpt_all(double ptmin,double ptmax);
    void compute_integrated_vn_all(double ptmin,double ptmax);

    // =============================
    // access
    // =============================
    event* get_event(int i) const noexcept { return event_arena[i]; }
    std::vector<event*> get_event_arena(){return event_arena;};
    int get_total_events() const noexcept { return total_music_events; }
    inline int get_Nptbins() const noexcept { return music_pt_bins; }
    inline double get_pt_val_of_bin(int ii){return ptval[ii];}

    void initial( const std::vector<std::string>& paths);
    
private:

    std::vector<std::string> music_output_paths;

    int max_Nevents;
    int yflag;

    double rapmin,rapmax,ptmin,ptmax;

    int total_music_events = 0;
    int music_pt_bins = 0;

    std::vector<double> ptval;
    std::vector<event*> event_arena;

    // species list (NO charged PID=0 here)
    const std::vector<int> PIDLIST =
        {211,-211,321,-321,2212,-2212};

    // integration helpers
    double integrate_spectrum(
        const std::vector<double>& _pt,
        const std::vector<double>& f,
        double ptmin,double ptmax);
    double integrate_spectrum_weighted(
        const std::vector<double>& _pt,
        const std::vector<double>& f,
        const std::vector<double>& w,
        double ptmin,double ptmax);


};
