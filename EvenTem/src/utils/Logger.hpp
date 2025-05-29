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

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <fstream>
#include <iostream>
#include <mutex>
#include <string>

class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    void log(const std::string& message) {
        std::lock_guard<std::mutex> guard(log_mutex);
        log_file << message << std::endl;
    }

private:
    Logger() : log_file("debug.log", std::ios_base::app) {
        if (!log_file.is_open()) {
            std::cerr << "Failed to open log file" << std::endl;
        }
    }

    ~Logger() {
        if (log_file.is_open()) {
            log_file.close();
        }
    }

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    std::ofstream log_file;
    std::mutex log_mutex;
};

class TXT_Logger {
    public:
        static TXT_Logger& getInstance() {
            static TXT_Logger instance;
            return instance;
        }
    
        void log(const std::string& message) {
            std::lock_guard<std::mutex> guard(log_mutex);
            log_file << message << std::endl;
        }
    
    private:
        TXT_Logger() : log_file("output.txt", std::ios_base::app) {
            if (!log_file.is_open()) {
                std::cerr << "Failed to open txt output file" << std::endl;
            }
        }
    
        ~TXT_Logger() {
            if (log_file.is_open()) {
                log_file.close();
            }
        }
    
        TXT_Logger(const TXT_Logger&) = delete;
        TXT_Logger& operator=(const TXT_Logger&) = delete;
    
        std::ofstream log_file;
        std::mutex log_mutex;
    };

#endif // LOGGER_HPP