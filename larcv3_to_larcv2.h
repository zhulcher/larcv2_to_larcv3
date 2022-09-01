#ifndef __larcv3_to_larcv2_H
#define __larcv3_to_larcv2_H


// Larcv2 includes:
#include "larcv3/core/dataformat/IOManager.h"

// Larcv3 includes:
#include "larcv/core/DataFormat/IOManager.h"


class larcv3_to_larcv2
{
public:
    // larcv3_to_larcv2();
    // ~larcv3_to_larcv2();
    
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
    larcv3_to_larcv2* larcv3_to_larcv2_new() {return new larcv3_to_larcv2();}
    void larcv3_to_larcv2_add_in_file (larcv3_to_larcv2* larcv3_to_larcv2, const char * s) {larcv3_to_larcv2->add_in_file(s);}
    void larcv3_to_larcv2_set_out_file(larcv3_to_larcv2* larcv3_to_larcv2, const char * s) {larcv3_to_larcv2->set_out_file(s);}
    void larcv3_to_larcv2_initialize(larcv3_to_larcv2* larcv3_to_larcv2) {larcv3_to_larcv2->initialize();}
    void larcv3_to_larcv2_convert(larcv3_to_larcv2* larcv3_to_larcv2, int n_events, int n_skip) {larcv3_to_larcv2->convert(n_events, n_skip);}
}

#endif
