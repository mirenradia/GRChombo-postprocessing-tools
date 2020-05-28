/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SMALLDATAIOREADER_HPP
#define SMALLDATAIOREADER_HPP

#include <fstream>
#include <string>
#include <vector>

// A class to read files written using SmallDataIO.

class SmallDataIOReader
{
  public:
    // A struct for information about the structure of a SmallDataIO file
    struct file_structure_t
    {
        int num_blocks; // a block is separated by 2 blank lines
        std::vector<unsigned long long int>
            block_starts; // position offsets from the beginning of the file
        std::vector<int> num_data_rows;      // the number of data rows in each
                                             // block
        std::vector<int> num_header_rows;    // the number of header rows in
                                             // each block
        std::vector<int> num_coords_columns; // number of coord columns in
                                             // each block - assume constant
                                             // in each block
        int coords_width;                    // assume the same throughout file
        std::vector<int> num_data_columns;   // number of data columns in each
                                             // block

        int data_width;               // assume the same throughout file
        std::vector<int> num_columns; // sum of coord and data columns in each
                                      // block
    };

  protected:
    std::string m_filename;
    std::ifstream m_file;
    file_structure_t m_file_structure;
    bool m_structure_defined;

  public:
    // Constructor
    SmallDataIOReader();

    // Destructor
    ~SmallDataIOReader();

    // Opens the file and sets m_filename. Note that this does not determine the
    // file structure
    void open(const std::string &a_filename);

    // Closes the file
    void close();

    // Parses the file and determines its structure
    void determine_file_structure();

    // Set structure if known already (e.g. same as another file already
    // determined)
    void set_file_structure(const file_structure_t &a_file_structure);

    // File struture getter
    const file_structure_t &get_file_structure() const;

    // Get a column from a block (either coord or data)
    std::vector<double> get_column(int a_column, int a_block = 0);

    // Get a data column from a block
    std::vector<double> get_data_column(int a_data_column, int a_block = 0);

    // Get same column from all blocks
    std::vector<std::vector<double>> get_column_from_all_blocks(int a_column);

    // Get same data column from all blocks
    std::vector<std::vector<double>>
    get_data_column_from_all_blocks(int a_data_column);
};

#endif /* SMALLDATAIOREADER_HPP */
