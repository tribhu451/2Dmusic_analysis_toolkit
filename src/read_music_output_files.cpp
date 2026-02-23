#include "read_music_output_files.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <gsl/gsl_integration.h>

// =================================================
// constructor
// =================================================
read_music_output_files::read_music_output_files(
    const std::vector<std::string>& paths,
    int max_nevents,
    int ayflag,
    double arapmin,double arapmax,
    double aptmin,double aptmax
)
: music_output_paths(paths),
  max_Nevents(max_nevents),
  yflag(ayflag),
  rapmin(arapmin),
  rapmax(arapmax),
  ptmin(aptmin),
  ptmax(aptmax)
{
    std::ifstream file;
    std::string line;

    int temp_events=0;

     std::cout << "reading started" << std::endl ; 
    // ---------- first pass ----------
    for(const auto& path:music_output_paths){
        for(int ie=0;ie<max_Nevents;ie++){
      

            std::ostringstream fname;
            fname<<path<<"/outputs_"
                 <<std::setw(4)<<std::setfill('0')<<ie
                 <<"/Fvnpt-211_y_-0.5_0.5.dat";

            file.open(fname.str());
            if(!file.is_open()) continue;

            std::cout << "event = " << ie << "  " << fname.str() << std::endl ; 
            temp_events++;

            if(temp_events==1){
                std::getline(file,line);
                double pt,dummy;
                while(std::getline(file,line)){
                    std::istringstream iss(line);
                    if(!(iss>>pt>>dummy)) break;
                    ptval.push_back(pt);
                }
            }
            file.close();
        }
    }

    total_music_events=temp_events;
    music_pt_bins=ptval.size();

    std::cout << "Nptbins = " << music_pt_bins << std::endl ; 
    std::cout << "Nevents = " << total_music_events << std::endl ; 
    if(total_music_events==0){
        std::cerr<<"No events found\n";
        std::exit(EXIT_FAILURE);
    }

    for(int i=0;i<total_music_events;i++)
        event_arena.push_back(new event(music_pt_bins));

    // ---------- second pass ----------
    int iev=0;

    for(const auto& path:music_output_paths){
        for(int ie=0;ie<max_Nevents;ie++){

            bool found=false;

            for(int PID:PIDLIST){

                std::ostringstream fname;
                fname<<path<<"/outputs_"
                     <<std::setw(4)<<std::setfill('0')<<ie
                     <<"/Fvnpt-"<<PID<<"_y_-0.5_0.5.dat";

                file.open(fname.str());
                if(!file.is_open()) continue;

                found=true;
                std::getline(file,line);

                double pt,dn,v1c,v1s,v2c,v2s,v3c,v3s,v4c,v4s,dum;
                int ipt=0;
                event* ev=event_arena[iev];

                while(std::getline(file,line)&&ipt<music_pt_bins){

                    std::istringstream iss(line);
                    iss>>pt>>dn
                       >>v1c>>v1s>>v2c>>v2s>>v3c>>v3s>>v4c>>v4s>>dum;

                    ev->set_differential_vn(PID,0,0,ipt,dn);

                    ev->set_differential_vn(PID,1,0,ipt,v1c);
                    ev->set_differential_vn(PID,1,1,ipt,v1s);
                    ev->set_differential_vn(PID,2,0,ipt,v2c);
                    ev->set_differential_vn(PID,2,1,ipt,v2s);
                    ev->set_differential_vn(PID,3,0,ipt,v3c);
                    ev->set_differential_vn(PID,3,1,ipt,v3s);
                    ev->set_differential_vn(PID,4,0,ipt,v4c);
                    ev->set_differential_vn(PID,4,1,ipt,v4s);

                    ipt++;
                }
                file.close();
            }

            if(found) iev++;
        }
    }
}

// =================================================
// destructor
// =================================================
read_music_output_files::~read_music_output_files(){
    for(auto e:event_arena) delete e;
}


// =================================================
// integration helpers
// =================================================
double read_music_output_files::integrate_spectrum(
        const std::vector<double>& _pt,
        const std::vector<double>& f,
        double minpt,double maxpt)
{
 
    int npt = _pt.size() ; 
    double dndpt[npt];
    double pt[npt];
    for (int ipt = 0; ipt < npt; ipt++) {
        pt[ipt]    = _pt[ipt];
        dndpt[ipt] = f[ipt] * 2 * M_PI * pt[ipt] ;  // dN/dpt
    }
    gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
    gsl_spline *numspline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (numspline, pt , dndpt , npt);
    
    double num = gsl_spline_eval_integ(numspline, minpt, maxpt, numacc);

    gsl_spline_free (numspline);
    gsl_interp_accel_free (numacc);
    
    return num;
}

double read_music_output_files::integrate_spectrum_weighted(
        const std::vector<double>& _pt,
        const std::vector<double>& f,
        const std::vector<double>& w,
        double minpt,double maxpt)
{
    int npt = _pt.size() ; 
    double dndpt[npt];
    double pt[npt];
    for (int ipt = 0; ipt < npt; ipt++) {
        pt[ipt]    = _pt[ipt];
        dndpt[ipt] = w[ipt] * f[ipt] * 2 * M_PI * pt[ipt] ;  //  dN/dpt * w(pt)
    }
    gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
    gsl_spline *numspline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (numspline, pt , dndpt , npt);
    
    double num = gsl_spline_eval_integ(numspline, minpt, maxpt, numacc);

    gsl_spline_free (numspline);
    gsl_interp_accel_free (numacc);
    
    return num;
}



// =================================================
// charged differential vn
// =================================================
void read_music_output_files::compute_differential_vn_charged_hadron(){

    for(auto ev:event_arena){
        for(int ipt=0; ipt<music_pt_bins; ipt++){

            double total=0;
            double sum[5][2]={{0}};

            for(int PID:PIDLIST){

                double yield=
                    ev->get_pt_differential_vn(PID,0,0,ipt);

                total+=yield;

                for(int h=1;h<5;h++)
                for(int ri=0;ri<2;ri++)
                    sum[h][ri]+=yield*
                        ev->get_pt_differential_vn(PID,h,ri,ipt);
            }

            if(total<=0) continue;

            for(int h=1;h<5;h++)
            for(int ri=0;ri<2;ri++)
                ev->set_differential_vn(0,h,ri,ipt,
                    sum[h][ri]/total);

            ev->set_differential_vn(0,0,0,ipt,total);
        }
    }
}


// =================================================
// calculators
// =================================================


// =================================================
// species mean pt
// =================================================
double read_music_output_files::calc_meanpt(
        event* ev,int PID,double ptmin,double ptmax)
{
    std::vector<double> spec(music_pt_bins);

    for(int i=0;i<music_pt_bins;i++){
        spec[i]=ev->get_pt_differential_vn(PID,0,0,i);
    }

    double numer=integrate_spectrum_weighted(ptval,spec,ptval,ptmin,ptmax);
    double denom=integrate_spectrum(ptval,spec,ptmin,ptmax);
    if(denom<=0) return 1e-20;

    return numer/denom;
}

double read_music_output_files::calc_integrated_vn(
        event* ev,int PID,int h,int ri,
        double ptmin,double ptmax)
{
    std::vector<double> spec(music_pt_bins);
    std::vector<double> vn(music_pt_bins);

    for(int i=0;i<music_pt_bins;i++){
        spec[i]=ev->get_pt_differential_vn(PID,0,0,i);
        vn[i]=ev->get_pt_differential_vn(PID,h,ri,i);
    }

    double numer=integrate_spectrum_weighted(ptval,spec,vn,ptmin,ptmax);
    double denom=integrate_spectrum(ptval,spec,ptmin,ptmax);
    if(denom<=0) return 0;

    if(h==0) 
      return denom; 
    else 
      return numer/denom;
}


// =================================================
//  setters
// =================================================
void read_music_output_files::compute_meanpt_all(
        double ptmin,double ptmax)
{
    for(auto ev:event_arena){
            ev->set_mean_pt(0,
                calc_meanpt(ev,0,ptmin,ptmax));
        for(int PID:PIDLIST){
            ev->set_mean_pt(PID,
                calc_meanpt(ev,PID,ptmin,ptmax));
        }};
}



void read_music_output_files::compute_integrated_vn_all(
        double ptmin,double ptmax)
{
    for(auto ev:event_arena){
        for(int PID:PIDLIST){
            for(int h=0;h<5;h++){
                for(int ri=0;ri<2;ri++){
                    ev->set_integrated_vn(PID,h,ri,
                        calc_integrated_vn(ev,PID,h,ri,ptmin,ptmax));
         }}}};

    for(auto ev:event_arena){
            for(int h=0;h<5;h++){
                for(int ri=0;ri<2;ri++){
                    ev->set_integrated_vn(0,h,ri,
                        calc_integrated_vn(ev,0,h,ri,ptmin,ptmax));
         }}};
}




