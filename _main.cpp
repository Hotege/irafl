
#if defined(_WIN32) && defined(_WINDLL)
#include <Windows.h>

BOOL WINAPI DllMain(HINSTANCE hModule, DWORD dwReason, LPVOID lpParam)
{
    UNREFERENCED_PARAMETER(hModule);
    UNREFERENCED_PARAMETER(dwReason);
    UNREFERENCED_PARAMETER(lpParam);
    switch (dwReason)
    {
    case DLL_PROCESS_ATTACH:
        break;
    case DLL_PROCESS_DETACH:
        break;
    case DLL_THREAD_ATTACH:
        break;
    case DLL_THREAD_DETACH:
        break;
    default:
        break;
    }
    return TRUE;
}

#endif