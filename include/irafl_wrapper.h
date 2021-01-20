#if !defined(_IRAFL_WRAPPER_H_)
#define _IRAFL_WRAPPER_H_

#include <stdint.h>

#if defined(__cplusplus)
extern "C" {
#endif

extern void IRAFL_Initialize();
extern void IRAFL_Terminate();

extern const char* IRAFL_Version();
extern const char* IRAFL_Info();

extern void IRAFL_Execute(unsigned char* dst, const unsigned char* src, long w, long h, uint32_t id, const void* params);
extern void IRAFL_Free(unsigned char* ptr);

typedef struct IRAFL_VP
{
    double r;
    uint32_t k;
} IRAFL_VP;

typedef struct IRAFL_LUT3DIMAGE
{
    unsigned char* data;
    long w, h, size;
} IRAFL_LUT3DIMAGE;

#if defined(__cplusplus)
};
#endif

#endif
