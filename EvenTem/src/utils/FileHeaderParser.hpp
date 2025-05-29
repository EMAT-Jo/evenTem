#include "FileConnector.h"

int parse_npy_header(std::string &file_path)
{

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
    
    std::vector<int> shape;
    size_t pos = 0;
    while ((pos = shape_str.find(',')) != std::string::npos) {
        shape.push_back(std::stoi(shape_str.substr(0, pos)));
        shape_str.erase(0, pos + 1);
    }
    if (!shape_str.empty()) {
        shape.push_back(std::stoi(shape_str));  // Add the last dimension
    }

    // Parse dtype
    std::string dtype;
    int bitdepth;
    dtype = header.substr(find_descr + 10, header.find("'", find_descr + 10) - (find_descr + 10));
    if (dtype == "|u1") {
        bitdepth = 8;
    } else if (dtype == "<u2") {
        bitdepth = 16;
    } else {
        throw std::runtime_error("Unsupported dtype");
    }

    file.close_file()

    return bitdepth;
};

void parse_mib_header()
{

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
            ds_merlin = stoi(head[4]) * stoi(head[5]);
            this->dtype = head[6];
            std::cout << "dtype: " << this->dtype << std::endl;
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
    if (this->dtype == "U08")
            {
                return 8;
            }
            else if (this->dtype == "U16")
            {
                return 16;
            }

}