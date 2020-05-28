/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIOReader.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    SmallDataIOReader file_reader;

    // open file
    std::cout << "Opening File\n";
    file_reader.open("test.dat");

    // determine file structure
    std::cout << "Determining File structure\n";
    file_reader.determine_file_structure();

    // Get file structure
    auto file_structure = file_reader.get_file_structure();
    std::cout << "num_blocks = " << file_structure.num_blocks << "\n";
    std::cout << "block_starts = ";
    for (auto start : file_structure.block_starts)
    {
        std::cout << start << " ";
    }
    std::cout << "\nnum_header_rows = ";
    for (auto num : file_structure.num_header_rows)
    {
        std::cout << num << " ";
    }
    std::cout << "\nnum_data_rows = ";
    for (auto num : file_structure.num_data_rows)
    {
        std::cout << num << " ";
    }
    std::cout << "\nnum_coords_columns = ";
    for (auto num : file_structure.num_coords_columns)
    {
        std::cout << num << " ";
    }
    std::cout << "\ncoords_width = " << file_structure.coords_width;
    std::cout << "\nnum_data_columns = ";
    for (auto num : file_structure.num_data_columns)
    {
        std::cout << num << " ";
    }
    std::cout << "\ndata_width = " << file_structure.data_width << std::endl;

    file_reader.get_column(2, 1);
}
