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
    void convert(int n_events = -1);

    void add_in_file(std::string filename);

    void set_out_file(std::string filename);

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

#endif
