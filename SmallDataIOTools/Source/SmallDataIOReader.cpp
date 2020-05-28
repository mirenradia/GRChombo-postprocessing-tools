/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIOReader.hpp"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

// A class to read files written using SmallDataIO.

// Constructor
SmallDataIOReader::SmallDataIOReader() : m_structure_defined(false) {}

// Destructor
SmallDataIOReader::~SmallDataIOReader()
{
    if (m_file.is_open())
    {
        close();
    }
}

// Opens the file and sets m_filename. Note that this does not determine the
// file structure
void SmallDataIOReader::open(const std::string &a_filename)
{
    m_filename = a_filename;
    m_structure_defined = false;
    m_file.open(m_filename);

    // check file opening successful
    if (!m_file)
    {
        std::cerr << "Error in opening " << m_filename << ". Exiting..."
                  << std::endl;
        exit(1);
    }
}

// Closes the file
void SmallDataIOReader::close()
{
    if (m_file.is_open())
    {
        m_file.close();
    }
    m_filename = "";
    m_structure_defined = false;
}

// Parses the file and determines its structure
void SmallDataIOReader::determine_file_structure()
{
    // first check file is open
    assert(m_file.is_open());

    // move file stream position to start of file
    m_file.clear();
    m_file.seekg(0, std::ios::beg);

    // go through each line and determine structure
    std::string line;
    int consecutive_empty_line_count = 2;
    int block_counter = 0; // assume we always have one block
    int current_position = m_file.tellg();
    int header_row_counter = 0;
    int data_row_counter = 0;
    while (std::getline(m_file, line))
    {
        if (!line.empty())
        {
            if (consecutive_empty_line_count == 2)
            {
                // start of new block
                ++block_counter;
                m_file_structure.block_starts.push_back(current_position);
            }
            consecutive_empty_line_count = 0;

            // header rows start with '#'
            if (line.find("#") != std::string::npos)
                ++header_row_counter;
            else
            {
                if (data_row_counter++ == 0)
                {
                    // determine column structure from first data row in block
                    // get a vector of the widths of the columns including
                    // preceeding whitespace
                    std::vector<int> widths;
                    std::string::size_type start_whitespace = 0;
                    while (!(start_whitespace == std::string::npos))
                    {
                        std::string::size_type start_non_whitespace =
                            line.find_first_not_of(' ', start_whitespace);
                        std::string::size_type next_start_whitespace =
                            line.find_first_of(' ', start_non_whitespace);
                        int width;
                        if (next_start_whitespace == std::string::npos)
                        {
                            width = line.length() - start_whitespace;
                        }
                        else
                        {
                            width = next_start_whitespace - start_whitespace;
                        }
                        widths.push_back(width);
                        start_whitespace = next_start_whitespace;
                    }

                    if (block_counter == 1)
                    {
                        // first data row in file so get coord and data width
                        // from this. assume min width is coord width and max is
                        // data width
                        auto widths_minmax_it =
                            std::minmax_element(widths.begin(), widths.end());
                        m_file_structure.coords_width =
                            *(widths_minmax_it.first);
                        m_file_structure.data_width =
                            *(widths_minmax_it.second);
                    }
                    int num_coords_columns =
                        std::count(widths.begin(), widths.end(),
                                   m_file_structure.coords_width);
                    int num_data_columns =
                        std::count(widths.begin(), widths.end(),
                                   m_file_structure.data_width);
                    m_file_structure.num_coords_columns.push_back(
                        num_coords_columns);
                    m_file_structure.num_data_columns.push_back(
                        num_data_columns);
                    m_file_structure.num_columns.push_back(widths.size());
                }
            }
        }
        else
        {
            if (consecutive_empty_line_count++ == 0 &&
                (header_row_counter > 0 || data_row_counter > 0))
            {
                // end of previous block
                m_file_structure.num_header_rows.push_back(header_row_counter);
                m_file_structure.num_data_rows.push_back(data_row_counter);
                header_row_counter = 0;
                data_row_counter = 0;
            }
        }
        current_position = m_file.tellg();
    }

    m_file_structure.num_blocks = block_counter;

    m_structure_defined = true;
}

// Set structure if known already (e.g. same as another file already
// determined)
void SmallDataIOReader::set_file_structure(
    const SmallDataIOReader::file_structure_t &a_file_structure)
{
    m_file_structure = a_file_structure;
}

// File struture getter
const SmallDataIOReader::file_structure_t &
SmallDataIOReader::get_file_structure() const
{
    return m_file_structure;
}

// Get a column from a block (either coord or data)
std::vector<double> SmallDataIOReader::get_column(int a_column, int a_block)
{
    assert(m_structure_defined);
    assert(a_column < m_file_structure.num_columns[a_block]);
    std::vector<double> out(m_file_structure.num_data_rows[a_block]);

    // how many characters across is the start of the column
    int start_position, column_width;
    if (a_column < m_file_structure.num_coords_columns[a_block])
    {
        start_position = a_column * m_file_structure.coords_width;
        column_width = m_file_structure.coords_width;
    }
    else
    {
        start_position =
            m_file_structure.num_coords_columns[a_block] *
                m_file_structure.coords_width +
            (a_column - m_file_structure.num_coords_columns[a_block]) *
                m_file_structure.data_width;
        column_width = m_file_structure.data_width;
    }

    // move stream position to start of block
    m_file.clear();
    m_file.seekg(m_file_structure.block_starts[a_block], std::ios::beg);
    std::string line;

    // assume header rows are all at the top of the block so skip these
    for (int irow = 0; irow < m_file_structure.num_header_rows[a_block]; ++irow)
    {
        std::getline(m_file, line);
    }
    for (int irow = 0; irow < m_file_structure.num_data_rows[a_block]; ++irow)
    {
        std::getline(m_file, line);
        out[irow] = std::stod(line.substr(start_position, column_width));
    }
    return out;
}

// Get a data column from a block
std::vector<double> SmallDataIOReader::get_data_column(int a_data_column,
                                                       int a_block)
{
    assert(m_structure_defined);
    int column = m_file_structure.num_coords_columns[a_block] + a_data_column;
    return get_column(column, a_block);
}

// Get same column from all blocks
std::vector<std::vector<double>>
SmallDataIOReader::get_column_from_all_blocks(int a_column)
{
    assert(m_structure_defined);
    std::vector<std::vector<double>> out(m_file_structure.num_blocks);

    for (int iblock = 0; iblock < m_file_structure.num_blocks; ++iblock)
    {
        out[iblock] = std::move(get_column(a_column, iblock));
    }
    return out;
}

// Get same data column from all blocks
std::vector<std::vector<double>>
SmallDataIOReader::get_data_column_from_all_blocks(int a_data_column)
{
    assert(m_structure_defined);
    std::vector<std::vector<double>> out(m_file_structure.num_blocks);

    for (int iblock = 0; iblock < m_file_structure.num_blocks; ++iblock)
    {
        out[iblock] = std::move(get_data_column(a_data_column, iblock));
    }
    return out;
}
