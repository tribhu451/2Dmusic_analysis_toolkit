#include <iostream>
#include <vector>

#include "event.h"
#include "read_music_output_files.h"
#include "observables.h"

int main(int argc,char** argv)
{


    int yflag = 1;

    double rapmin=-0.5;
    double rapmax= 0.5;

    double ptmin=0.2;
    double ptmax=3.0;

    int max_Nevents=999;

    if(argc<2){
        std::cout<<"Need at least one path\n";
        return -1;
    }

    std::vector<std::string> paths;
    for(int i=1;i<argc;i++)
        paths.push_back(argv[i]);

    read_music_output_files* rmof =
        new read_music_output_files(
            paths,max_Nevents,
            yflag,
            rapmin,rapmax,
            ptmin,ptmax);

    // build charged differential first
    rmof->compute_differential_vn_charged_hadron();

    // species
    rmof->compute_integrated_vn_all(ptmin,ptmax);
    rmof->compute_meanpt_all(ptmin,ptmax);
    //rmof-> initial(paths);
     
    observables* obj =
        new observables(rmof,yflag,
                        rapmin,rapmax,
                        ptmin,ptmax);

    obj->output_dndy_or_dndeta();
    
    
    obj->output_vo_ratio_proton_hpm();
    obj->output_ebe_meanpt_correlation_hpm_proton();
    obj-> output_Bozek_rn(2, 0, 2212); // harmonics, PID1, PID2
    obj-> output_Bozek_rn(3, 0, 2212);
    obj-> output_vn_meanpt_correlation_ratio_proton_charged_hadron(2);
    obj-> output_vn_meanpt_correlation_ratio_proton_charged_hadron(3);

    /*
    obj->output_pt_diff_multiparticle_vn_charged_hadrons(2);
    obj->output_pt_diff_multiparticle_vn_charged_hadrons(3);

    obj->output_pt_diff_multiparticle_vn_method2_charged_hadrons(2);
    obj->output_pt_diff_multiparticle_vn_method2_charged_hadrons(3);
    */
    return 0;
}



