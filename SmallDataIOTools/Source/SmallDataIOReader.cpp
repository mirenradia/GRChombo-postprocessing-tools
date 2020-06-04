/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIOReader.hpp"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <regex>

void SmallDataIOReader::file_structure_t::clear()
{
    num_blocks = 0;
    block_starts.clear();
    num_data_rows.clear();
    num_header_rows.clear();
    num_coords_columns.clear();
    coords_width = 0;
    num_data_columns.clear();
    data_width = 0;
    num_columns.clear();
}

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
    m_filename.clear();
    m_structure_defined = false;
    m_file_structure.clear();
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
    // Just in case the file ends without a line break:
    if (header_row_counter > 0 || data_row_counter > 0)
    {
        m_file_structure.num_header_rows.push_back(header_row_counter);
        m_file_structure.num_data_rows.push_back(data_row_counter);
        header_row_counter = 0;
        data_row_counter = 0;
    }

    m_file_structure.num_blocks = block_counter;

    assert(m_file_structure.num_data_rows.size() ==
               m_file_structure.num_blocks &&
           m_file_structure.num_header_rows.size() ==
               m_file_structure.num_blocks &&
           m_file_structure.num_coords_columns.size() ==
               m_file_structure.num_blocks &&
           m_file_structure.num_data_columns.size() ==
               m_file_structure.num_blocks &&
           m_file_structure.num_columns.size() == m_file_structure.num_blocks &&
           m_file_structure.block_starts.size() == m_file_structure.num_blocks);

    m_structure_defined = true;
}

// Set structure if known already (e.g. same as another file already
// determined)
void SmallDataIOReader::set_file_structure(
    const SmallDataIOReader::file_structure_t &a_file_structure)
{
    m_file_structure = a_file_structure;
    m_structure_defined = true;
}

// File struture getter
const SmallDataIOReader::file_structure_t &
SmallDataIOReader::get_file_structure() const
{
    return m_file_structure;
}

// Get an interval of columns (inclusive) from a block
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_columns(int a_min_column, int a_max_column, int a_block)
{
    assert(m_file.is_open());
    assert(m_structure_defined);
    assert(0 <= a_min_column <= a_max_column <
           m_file_structure.num_columns[a_block]);
    const int num_columns = a_max_column - a_min_column + 1;
    std::vector<column_t> out(num_columns);
    for (auto &column : out)
    {
        column.resize(m_file_structure.num_data_rows[a_block]);
    }

    // how many characters across is the start of each column
    std::vector<int> start_position(num_columns);
    std::vector<int> column_widths(num_columns);
    for (int ifcolumn = a_min_column; ifcolumn <= a_max_column; ++ifcolumn)
    {
        // ifcolumn is column index in file and icolumn is column index in
        // output
        int icolumn = ifcolumn - a_min_column;
        if (ifcolumn < m_file_structure.num_coords_columns[a_block])
        {
            start_position[icolumn] = ifcolumn * m_file_structure.coords_width;
            column_widths[icolumn] = m_file_structure.coords_width;
        }
        else
        {
            start_position[icolumn] =
                m_file_structure.num_coords_columns[a_block] *
                    m_file_structure.coords_width +
                (ifcolumn - m_file_structure.num_coords_columns[a_block]) *
                    m_file_structure.data_width;
            column_widths[icolumn] = m_file_structure.data_width;
        }
    }

    // move stream position to start of block
    m_file.clear();
    m_file.seekg(m_file_structure.block_starts[a_block], std::ios::beg);
    std::string line;

    // assume header rows are all at the top of the block so skip these
    for (int irow = 0; irow < m_file_structure.num_header_rows[a_block]; ++irow)
    {
        m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    for (int irow = 0; irow < m_file_structure.num_data_rows[a_block]; ++irow)
    {
        std::getline(m_file, line);
        for (int icolumn = 0; icolumn < num_columns; ++icolumn)
        {
            out[icolumn][irow] = std::stod(
                line.substr(start_position[icolumn], column_widths[icolumn]));
        }
    }
    return out;
}

// Get a data column from a block
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_all_data_columns(int a_block)
{
    assert(m_structure_defined);
    int min_data_column = m_file_structure.num_coords_columns[a_block];
    int max_data_column =
        min_data_column + m_file_structure.num_data_columns[a_block] - 1;
    return get_columns(min_data_column, max_data_column, a_block);
}

SmallDataIOReader::column_t SmallDataIOReader::get_column(int a_column,
                                                          int a_block)
{
    auto out_vect = std::move(get_columns(a_column, a_column, a_block));
    return out_vect[0];
}

// Returns a vector of numeric values from a header row
std::vector<double>
SmallDataIOReader::get_data_from_header(int a_header_row_number, int a_block)
{
    assert(m_file.is_open());
    assert(m_structure_defined);
    assert(a_header_row_number < m_file_structure.num_header_rows[a_block]);

    // move stream to start of block
    m_file.clear();
    m_file.seekg(m_file_structure.block_starts[a_block], std::ios::beg);
    std::string line;

    // assume header rows are at start of block
    for (int irow = 0; irow < a_header_row_number; ++irow)
    {
        // skip lines before desired row
        m_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // get desired header line
    std::getline(m_file, line);

    // find numbers in header using regex
    // I think this takes a long time to compile...
    std::regex number("[+-]?([0-9]*\\.)?[0-9]+");
    auto numbers_begin = std::sregex_iterator(line.begin(), line.end(), number);
    auto numbers_end = std::sregex_iterator();
    std::vector<double> out;
    for (std::sregex_iterator rit = numbers_begin; rit != numbers_end; ++rit)
    {
        // put matches in vector
        out.push_back(std::stod((*rit).str()));
    }

    return out;
}
