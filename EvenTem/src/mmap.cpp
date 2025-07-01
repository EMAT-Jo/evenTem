#include <string> 
#include <windows.h>
#include <iostream>
#include <cstdint>
#include <chrono>
#include <fstream>
#include <array>
#define min(a,b) (((a) < (b)) ? (a) : (b))


#ifdef __GNUC__
#define PACK(__Declaration__) __Declaration__ __attribute__((__packed__))
#endif

#ifdef _MSC_VER
#define PACK(__Declaration__) __pragma(pack(push, 1)) __Declaration__ __pragma(pack(pop))
#endif

PACK(struct EVENT
{
    uint32_t index;
    uint64_t toa;
    uint8_t overflow;
    uint8_t ftoa;
    uint16_t tot;
});

#define BUFFER_SIZE 512 * 1024 * 1024 // 512 MB
#define BUFFER_SIZE_EVENTS 65536 // Number of events that can be stored in the buffer


int main() {

    std::array<uint64_t,256*256> pacbed;
    for (size_t i = 0; i < pacbed.size(); ++i) {
        pacbed[i] = 0;
    }

    //---------------mmap-------------------
    
    std::string filename_str = "c:\\DATA-C\\STO_Christoph\\Mon_Feb_26_11_21_23_2024_STEM_200kV_1024x1024_5.0_us_5_scans.t3p";
    const char* filename = filename_str.c_str();

    HANDLE hFile = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    HANDLE hMap = CreateFileMappingA(hFile, NULL, PAGE_READONLY, 0, 0, NULL);

    size_t map_size = BUFFER_SIZE; 
    size_t offset = 0;
    LARGE_INTEGER size;
    GetFileSizeEx(hFile, &size);
    size_t file_size = size.QuadPart;
    std::cout << "File size: " << file_size << " bytes." << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    void* mapped_ptr = MapViewOfFile(hMap,FILE_MAP_READ,DWORD(offset >> 32),DWORD(offset & 0xFFFFFFFF),map_size);
    int total_num_events = 0;
    int num_events = 0;
    EVENT* events = static_cast<EVENT*>(mapped_ptr);
    offset += map_size;
    num_events += map_size / sizeof(EVENT);

    while (offset < file_size) 
    {
        if (file_size - offset < map_size) {
            map_size = file_size - offset;
            num_events = map_size / sizeof(EVENT);
        }

        mapped_ptr = MapViewOfFile(hMap,FILE_MAP_READ,DWORD(offset >> 32),DWORD(offset & 0xFFFFFFFF),map_size);

        if (mapped_ptr == NULL) {
            std::cerr << "Error mapping file: " << GetLastError() << std::endl;
            CloseHandle(hMap);
            CloseHandle(hFile);
            return 1;
        }

        events = static_cast<EVENT*>(mapped_ptr);
        for (int i = 0; i < num_events; ++i) {
            pacbed[events[i].index]++;
        }

        offset += map_size;
        total_num_events += map_size / sizeof(EVENT);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "mmap: "<< elapsed.count() << " seconds." << std::endl;
    std::cout << "Number of events: " << total_num_events << std::endl;

    // Cleanup
    UnmapViewOfFile(mapped_ptr);
    CloseHandle(hMap);
    CloseHandle(hFile);

    // for (size_t i = 0; i < pacbed.size(); ++i) {
    //     if (pacbed[i] > 0) {
    //         std::cout << "Index: " << i << ", Count: " << pacbed[i] << std::endl;
    //     }
    // }
    std::cout << pacbed[4242] << std::endl; // Example output for a specific index

    for (size_t i = 0; i < pacbed.size(); ++i) {
            pacbed[i] = 0;
        }

    //---------------compare to ifstream-------------------
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file with ifstream." << std::endl;
        CloseHandle(hMap);
        CloseHandle(hFile);
        return 1;
    }

    file.seekg(0, std::ios::end);
    size_t file_size_ifstream = file.tellg();
    file.seekg(0, std::ios::beg);
    std::cout << "File size: " << file_size_ifstream << " bytes." << std::endl;

    total_num_events = 0;

    std::array<EVENT, BUFFER_SIZE_EVENTS> *buffer = new std::array<EVENT, BUFFER_SIZE_EVENTS>;

    size_t bytes_to_read = sizeof(*buffer);
    auto start_time_ifstream = std::chrono::high_resolution_clock::now();

    size_t bytes_read = 0;

    while (bytes_read < file_size_ifstream) {
        bytes_to_read = min(file_size_ifstream - bytes_read, sizeof(*buffer));
        file.read(reinterpret_cast<char*>(buffer->data()), bytes_to_read);
        bytes_read += bytes_to_read;
        total_num_events += bytes_to_read / sizeof(EVENT);

        for (const auto& event : *buffer) {
            pacbed[event.index]++;
        }
    }
    
    auto end_time_ifstream = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ifstream = end_time_ifstream - start_time_ifstream;
    std::cout << "ifstream: " << elapsed_ifstream.count() << " seconds." << std::endl;
    std::cout << "Number of events: " << total_num_events << std::endl;

    std::cout << pacbed[4242] << std::endl;

    return 0;
}