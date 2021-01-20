#include <irafl_wrapper.h>
#include <irafl.h>

#if defined(__cplusplus)
extern "C" {
#endif

void IRAFL_Initialize()
{
    IRAFLH->Initialize();
}

void IRAFL_Terminate()
{
    IRAFLH->Terminate();
}

const char* IRAFL_Version()
{
    return IRAFL::Version();
}

const char* IRAFL_Info()
{
    return IRAFL::Info();
}

void IRAFL_Execute(unsigned char* dst, const unsigned char* src, long w, long h, uint32_t id, const void* params)
{
    IRAFL::Execute(dst, src, w, h, id, params);
}

void IRAFL_Free(unsigned char* ptr)
{
    delete[]ptr;
    ptr = nullptr;
}

#if defined(__cplusplus)
};
#endif
