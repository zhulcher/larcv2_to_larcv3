

#include <iostream>

#include "larcv2_to_larcv3.h"

int main()
{

    // Create a converter object
    larcv2_to_larcv3 converter; 

    // For development, hardcode input and output files:
    std::string input_file = "./dlprod_unified_00.root";
    std::string output_file = "./dlprod_unified_00.h5";

    converter.add_in_file(input_file);
    converter.set_out_file(output_file);

    converter.initialize(); 

    std::cout << "Converting!" << std::endl;

    converter.convert(200);

    return 0;
}