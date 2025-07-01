/* Copyright (C) 2025 Thomas Friedrich, Chu-Ping Yu, Arno Annys
 * University of Antwerp - All Rights Reserved. 
 * You may use, distribute and modify
 * this code under the terms of the GPL3 license.
 * You should have received a copy of the GPL3 license with
 * this file. If not, please visit: 
 * https://www.gnu.org/licenses/gpl-3.0.en.html
 * 
 * Authors: 
 *   Thomas Friedrich <>
 *   Chu-Ping Yu <>
 *   Arno Annys <arno.annys@uantwerpen.be>
 */

#ifndef FILE_CONNECTOR_MMAP_H
#define FILE_CONNECTOR_MMAP_H

#ifdef _WIN32
#include <windows.h>
#endif
#include <limits>
#include <string> 
#include <fstream>


class FileConnectorMmap
{
public:
    std::string path;
    std::uintmax_t file_size;
    std::uintmax_t pos;
    void open_file();
    void close_file();
    void read_data(void*& mapped_ptr, size_t map_size);
    FileConnectorMmap();

private:
    #ifdef _WIN32
    HANDLE hFile;
    HANDLE hMap;
    #endif


    void reset_file();
};
#endif // FILE_CONNECTOR_MMAP_H