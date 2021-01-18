#if !defined(_FUNCTIONS_GAUSSBLUR_HPP_)
#define _FUNCTIONS_GAUSSBLUR_HPP_

#include <irafldef.h>
#include "functions.hpp"
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>

namespace IRAFLFunc
{
    namespace name
    {
        static std::string GaussBlurName = "gaussblur";
    }
    namespace gaussblur
    {
        void IRAFL_LOCAL gaussblur(unsigned char*& dst, const unsigned char* src, const long& w, const long& h, const void* params)
        {
            memset(dst, 0, w * h * 4);
            long r = *(long*)params;
            long d = 2 * r + 1;
            std::vector<double> kernel(d, 0);
            double sigma = r / 3.0;
            double sigma2 = 2 * sigma * sigma;
            double sigmap = sigma * sqrt(2 * PI);
            for (long n = 0, i = -r; i <= r; ++i, ++n)
                kernel[n] = exp(-((double)i * i) / sigma2) / sigmap;
            double sum = 0;
            for (auto p : kernel)
                sum += p;
            for (size_t i = 0; i < kernel.size(); i++)
                kernel[i] /= sum;
            double* cr = new double[w * h * 4];
            memset(cr, 0, sizeof(double) * w * h * 4);
            for (long y = 0; y < h; y++)
                for (long x = 0; x < w; x++)
                    for (long i = -r; i <= r; i++)
                    {
                        long ii = x + i;
                        if (x + i < 0)
                            ii = 0;
                        else if (x + i >= w)
                            ii = w - 1;
                        cr[(y * w + x) * 4 + 0] += kernel[i + r] * src[(y * w + ii) * 4 + 0];
                        cr[(y * w + x) * 4 + 1] += kernel[i + r] * src[(y * w + ii) * 4 + 1];
                        cr[(y * w + x) * 4 + 2] += kernel[i + r] * src[(y * w + ii) * 4 + 2];
                    }
            for (long x = 0; x < w; x++)
                for (long y = 0; y < h; y++)
                {
                    double B = 0, G = 0, R = 0;
                    for (long j = -r; j <= r; j++)
                    {
                        long jj = y + j;
                        if (y + j < 0)
                            jj = 0;
                        else if (y + j >= h)
                            jj = h - 1;
                        B += kernel[j + r] * cr[(jj * w + x) * 4 + 0];
                        G += kernel[j + r] * cr[(jj * w + x) * 4 + 1];
                        R += kernel[j + r] * cr[(jj * w + x) * 4 + 2];
                    }
                    dst[(y * w + x) * 4 + 0] = (unsigned char)(B);
                    dst[(y * w + x) * 4 + 1] = (unsigned char)(G);
                    dst[(y * w + x) * 4 + 2] = (unsigned char)(R);
                    dst[(y * w + x) * 4 + 3] = 255;
                }
            delete[]cr;
            cr = nullptr;
            return;
        }
    }
}

#endif
