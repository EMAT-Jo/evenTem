
#ifndef FILEHEADERPARSER_HPP
#define FILEHEADERPARSER_HPP

#ifdef _WIN32
#include <io.h>
#pragma warning(disable : 4005 4333 34)
#else
#include <unistd.h>
#endif

#include <iostream>
#include <limits>
#include <string>
#include <fstream>
#include <filesystem> 
#include <array>
#include <sstream>
#include "FileConnector.h"

inline std::array<int,5> parse_npy_header(std::string &file_path)
{
    std::array<int, 5> header_info = {0, 0, 0, 0, 0}; // Initialize with zeros

    FileConnector file;
    file.path = file_path;
    file.open_file();

    // Read magic string
    char magic[6];
    file.read_data(magic, 6);
    if (std::memcmp(magic, "\x93NUMPY", 6) != 0) {
        throw std::runtime_error("Not a valid .npy file");
    }

    // Read version number
    uint8_t major_version, minor_version;
    file.read_data(reinterpret_cast<char*>(&major_version), 1);
    file.read_data(reinterpret_cast<char*>(&minor_version), 1);

    int16_t header_len;
    if (major_version == 1){
        file.read_data(reinterpret_cast<char*>(&header_len), 2);
        header_len = static_cast<int>(header_len);        
    }
    else{
        throw std::runtime_error("Unsupported version number of npy file, should be version 1 for simple structures");
    }

    // Read header
    std::string header(header_len, ' ');
    file.read_data(&header[0], header_len);

    // Extract shape and dtype
    auto find_shape = header.find("'shape': (");
    auto find_descr = header.find("'descr': '");
    auto find_order = header.find("'fortran_order': ");

    if (find_shape == std::string::npos || find_descr == std::string::npos) {
        throw std::runtime_error("Header parsing failed");
    }

    // Parse shape
    size_t start = find_shape + 10;
    size_t end = header.find(")", start);
    std::string shape_str = header.substr(start, end - start);
    
    size_t pos = 0;
    int count = 0;
    while ((pos = shape_str.find(',')) != std::string::npos) {
        header_info[count] = std::stoi(shape_str.substr(0, pos));
        count++;
        shape_str.erase(0, pos + 1);
    }
    if (!shape_str.empty()) {
        header_info[count] = std::stoi(shape_str);
    }

    // Parse dtype
    std::string dtype;
    dtype = header.substr(find_descr + 10, header.find("'", find_descr + 10) - (find_descr + 10));
    if (dtype == "|u1") {
        header_info[4] = 8;
    } else if (dtype == "<u2") {
        header_info[4] = 16;
    } else {
        throw std::runtime_error("Unsupported dtype");
    }

    file.close_file();  

    return header_info;
};

inline std::array<int,5> parse_mib_header(std::string &file_path)
{
    std::array<int, 5> header_info = {0, 0, 0, 0, 0}; // Initialize with zeros

    FileConnector file;
    file.path = file_path;
    file.open_file();

    std::string dtype;

    std::string rcv;

    std::array<char, 384> head_buffer;
    std::array<std::string, 8> head;
    file.read_data(&head_buffer[0], head_buffer.size());
    rcv.assign(head_buffer.cbegin(), head_buffer.cend());
    size_t i = 0;
    head.fill("");
    std::stringstream ss(rcv);
    while (ss.good() && i <= 6)
    {
        std::getline(ss, head[i], ',');
        i++;
    }
    if (i >= 6)
    {
        try
        {
            header_info[2] = stoi(head[4]);
            header_info[3] = stoi(head[5]);
            dtype = head[6];
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
    }
    else
    {
        perror("Frame Header cannot be decoded!");
    }

    if (dtype == "U08")
    {
        header_info[4] = 8;
    }
    else if (dtype == "U16")
    {
        header_info[4] = 16;
    }
    else {
        throw std::runtime_error("Unsupported dtype");
    }

    file.close_file();  

    return header_info;
}

#endif // FILEHEADERPARSER_HPP
