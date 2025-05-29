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

#ifndef DTYPE_ELECTRON_HPP
#define DTYPE_ELECTRON_HPP

#ifdef __GNUC__
#define PACK(__Declaration__) __Declaration__ __attribute__((__packed__))
#endif

#ifdef _MSC_VER
#define PACK(__Declaration__) __pragma(pack(push, 1)) __Declaration__ __pragma(pack(pop))
#endif

#include <stdint.h>

#pragma pack(push, 1)
struct dtype_Electron
{
    uint16_t kx;
    uint16_t ky;
    uint16_t rx;
    uint16_t ry;
    uint16_t id_image;
};
#pragma pack(pop)

#endif // DTYPE_ELECTRON_HPP