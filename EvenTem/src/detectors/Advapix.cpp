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

#include "Advapix.hpp"

#ifdef PIXET_ENABLED
inline void OnTpx3Data(intptr_t eventData, intptr_t userData)
{
    auto* cbData = reinterpret_cast<CallbackData_t*>(userData);
    int deviceIndex = cbData->deviceIndex;
    int n_buffer = ADVAPIX_ADDITIONAL::N_BUFFER;

    if (!(*cbData->p_n_buffer_filled < (n_buffer + *cbData->p_n_buffer_processed)))
    {
        std::cout << "Error: Buffer is full" << std::endl;
        return;
    }

    unsigned pixelCount = 0;
    int rc = pxcGetMeasuredTpx3PixelsCount(deviceIndex, &pixelCount);
    printf("Counts: %u\n", pixelCount);

    if (rc) return;

    std::shared_ptr<Tpx3Pixel[]> Tpx3Pixels(new Tpx3Pixel[pixelCount], std::default_delete<Tpx3Pixel[]>());
    rc = pxcGetMeasuredTpx3Pixels(deviceIndex, Tpx3Pixels.get(), pixelCount);

    if (rc) return;

    (*cbData->p_ragged_buffer)[*cbData->p_n_buffer_filled % n_buffer] = Tpx3Pixels;
    (*cbData->p_ragged_buffer_sizes)[*cbData->p_n_buffer_filled % n_buffer] = pixelCount;
    (*cbData->p_n_buffer_filled)++;
    return;
};
#endif
