#ifndef __LARCV2_TO_LARCV3_H
#define __LARCV2_TO_LARCV3_H


// Larcv2 includes:
#include "larcv3/core/dataformat/IOManager.h"

// Larcv3 includes:
#include "larcv/core/DataFormat/IOManager.h"


class larcv2_to_larcv3
{
public:
    // larcv2_to_larcv3();
    // ~larcv2_to_larcv3();
    
    void initialize();
    void convert(int n_events = -1, int n_skip = 0);

    void add_in_file(const char * filename);

    void set_out_file(const char * filename);

private:

    void convert_event(size_t entry);

    // Conversion functions:
    void convert_image2d(std::string producer);

    void convert_particle(std::string producer);
    void convert_sparse2d(std::string producer);
    void convert_sparse3d(std::string producer);
    void convert_cluster2d(std::string producer);
    void convert_cluster3d(std::string producer);

    std::vector<std::string> input_files;
    std::string output_file;


    larcv::IOManager larcv2_manager;
    larcv3::IOManager larcv3_manager;

};

// Define C functions for the C++ class - as ctypes can only talk to C...

extern "C"
{
    larcv2_to_larcv3* larcv2_to_larcv3_new() {return new larcv2_to_larcv3();}
    void larcv2_to_larcv3_add_in_file (larcv2_to_larcv3* larcv2_to_larcv3, const char * s) {larcv2_to_larcv3->add_in_file(s);}
    void larcv2_to_larcv3_set_out_file(larcv2_to_larcv3* larcv2_to_larcv3, const char * s) {larcv2_to_larcv3->set_out_file(s);}
    void larcv2_to_larcv3_initialize(larcv2_to_larcv3* larcv2_to_larcv3) {larcv2_to_larcv3->initialize();}
    void larcv2_to_larcv3_convert(larcv2_to_larcv3* larcv2_to_larcv3, int n_events, int n_skip) {larcv2_to_larcv3->convert(n_events, n_skip);}
}

#endif
