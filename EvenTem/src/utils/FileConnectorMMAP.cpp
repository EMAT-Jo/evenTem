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

#include "FileConnectorMmap.h"
#include <iostream>

void FileConnectorMmap::open_file()
{
    const char* filename = path.c_str();
    #ifdef _WIN32
        hFile = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        hMap = CreateFileMappingA(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    #endif
}

void FileConnectorMmap::close_file()
{
    #ifdef _WIN32
    // UnmapViewOfFile(mapped_ptr);
        CloseHandle(hMap);
        CloseHandle(hFile);
    #endif
}

// Reading data stream from File
void FileConnectorMmap::read_data(void*& mapped_ptr, size_t map_size)
{
    #ifdef _WIN32
        mapped_ptr = MapViewOfFile(hMap,FILE_MAP_READ,DWORD(pos >> 32),DWORD(pos & 0xFFFFFFFF),map_size);
        pos += map_size;
    #endif
};

void FileConnectorMmap::reset_file()
{

}

FileConnectorMmap::FileConnectorMmap(): path(), file_size(0), pos(0) {};