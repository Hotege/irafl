#include <irafl.h>
#include <cstring>
#include <ctime>
#include <string>
#include <map>
#include "functions/functions.hpp"
#include "functions/gray.hpp"
#include "functions/gaussblur.hpp"
#include "functions/lut3d.hpp"
#include "functions/voronoi.hpp"

#define INSERT_FILTER(id, f, tag) inner->filters.insert(std::make_pair(id, InnerImpl::FilterUnit{ IRAFLFunc::f::f, IRAFLFunc::name::tag##Name }))

IRAFL* IRAFL::instance = nullptr;
const char* IRAFL::version = IRAFL_VERSION;
const char* IRAFL::info = IRAFL_INFO;

class IRAFL::InnerImpl
{
public:
    typedef void (*FilterFunc)(unsigned char*&, const unsigned char*, const long&, const long&, const void*);
    typedef struct _FilterUnit
    {
        FilterFunc func;
        std::string name;
    } FilterUnit;
    std::map<uint32_t, FilterUnit> filters;
    std::string none = "";
};

IRAFL::IRAFL()
{
    inner = new InnerImpl;
    INSERT_FILTER(IRAFL_FILTER_TYPE_GRAY,       gray,       Gray);
    INSERT_FILTER(IRAFL_FILTER_TYPE_GAUSSBLUR,  gaussblur,  GaussBlur);
    INSERT_FILTER(IRAFL_FILTER_TYPE_LUT3D,      lut3d,      Lut3D);
    INSERT_FILTER(IRAFL_FILTER_TYPE_VORONOI,    voronoi,    Voronoi);
}

IRAFL::~IRAFL()
{
    delete inner;
    inner = nullptr;
}

IRAFL* IRAFL::Handle()
{
    if (!instance)
        instance = new IRAFL;
    return instance;
}

void IRAFL::Initialize()
{
    unsigned int ct = time(nullptr);
    srand(IRAFLFunc::CRC32Get(&ct, sizeof(unsigned int)));
}

void IRAFL::Terminate()
{
    delete instance;
    instance = nullptr;
}

const char* IRAFL::Version()
{
    return version;
}

const char* IRAFL::Info()
{
    return info;
}

void IRAFL::Execute(unsigned char*& dst, const unsigned char* src, const long& w, const long& h, const uint32_t& id, const void* params)
{
    if (IRAFLH->inner->filters.count(id) != 0)
        IRAFLH->inner->filters[id].func(dst, src, w, h, params);
}

const std::string& IRAFL::GetFilterName(const uint32_t& id)
{
    if (IRAFLH->inner->filters.count(id) != 0)
        return IRAFLH->inner->filters[id].name;
    else
        return IRAFLH->inner->none;
}
