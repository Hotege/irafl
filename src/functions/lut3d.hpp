#if !defined(_FUNCTIONS_LUT3D_HPP_)
#define _FUNCTIONS_LUT3D_HPP_

#include <irafldef.h>
#include "functions.hpp"
#include <cstdint>
#include <string>

namespace IRAFLFunc
{
    namespace name
    {
        static std::string Lut3DName = "lut3d";
    }
    namespace lut3d
    {
        typedef struct _Lut3DImage
        {
            unsigned char* data;
            long w, h, size;
        } Lut3DImage;
        void IRAFL_LOCAL lut3d(unsigned char*& dst, const unsigned char* src, const long& w, const long& h, const void* params)
        {
            Lut3DImage* img = (Lut3DImage*)params;
            if (img->w != img->h)
                return;
            for (long y = 0; y < h; y++)
                for (long x = 0; x < w; x++)
                {
                    unsigned char B = src[(y * w + x) * 4 + 0];
                    unsigned char G = src[(y * w + x) * 4 + 1];
                    unsigned char R = src[(y * w + x) * 4 + 2];
                    double lR = ColorNormal[R], lG = ColorNormal[G], lB = ColorNormal[B];
                    long rf = floor(lR * (img->size - 1)), rl = ceil(lR * (img->size - 1));
                    long gf = floor(lG * (img->size - 1)), gl = ceil(lG * (img->size - 1));
                    long bf = floor(lB * (img->size - 1)), bl = ceil(lB * (img->size - 1));
                    double tR = lR * (img->size - 1) - rf;
                    double tG = lG * (img->size - 1) - gf;
                    double tB = lB * (img->size - 1) - bf;

                    long xFF = bf % (img->w / img->size) * img->size + rf;
                    long xFL = bf % (img->w / img->size) * img->size + rl;
                    long yFF = bf / (img->w / img->size) * img->size + gf;
                    long yFL = bf / (img->w / img->size) * img->size + gl;
                    double crF[3] = { 0 };
                    for (int k = 0; k < 3; k++)
                    {
                        unsigned char crFF = img->data[(yFF * img->w + xFF) * 4 + k];
                        unsigned char crLF = img->data[(yFF * img->w + xFL) * 4 + k];
                        unsigned char crFL = img->data[(yFL * img->w + xFF) * 4 + k];
                        unsigned char crLL = img->data[(yFL * img->w + xFL) * 4 + k];
                        crF[k] = (crFF * (1 - tR) + crLF * tR) * (1 - tG) + (crFL * (1 - tR) + crLL * tR) * tG;
                    }

                    long xLF = bl % (img->w / img->size) * img->size + rf;
                    long xLL = bl % (img->w / img->size) * img->size + rl;
                    long yLF = bl / (img->w / img->size) * img->size + gf;
                    long yLL = bl / (img->w / img->size) * img->size + gl;
                    double crL[3] = { 0 };
                    for (int k = 0; k < 3; k++)
                    {
                        unsigned char crFF = img->data[(yLF * img->w + xLF) * 4 + k];
                        unsigned char crLF = img->data[(yLF * img->w + xLL) * 4 + k];
                        unsigned char crFL = img->data[(yLL * img->w + xLF) * 4 + k];
                        unsigned char crLL = img->data[(yLL * img->w + xLL) * 4 + k];
                        crL[k] = (crFF * (1 - tR) + crLF * tR) * (1 - tG) + (crFL * (1 - tR) + crLL * tR) * tG;
                    }

                    dst[(y * w + x) * 4 + 0] = crF[0] * (1 - tB) + crL[0] * tB;
                    dst[(y * w + x) * 4 + 1] = crF[1] * (1 - tB) + crL[1] * tB;
                    dst[(y * w + x) * 4 + 2] = crF[2] * (1 - tB) + crL[2] * tB;
                    dst[(y * w + x) * 4 + 3] = src[(y * w + x) * 4 + 3];
                }
        }
    }
}

#endif
