#if !defined(_IRAFL_DEF_H_)
#define _IRAFL_DEF_H_

#if defined(_WIN32) || defined(__CYGWIN__)
#if defined(_WINDLL) || defined(BUILDING_DLL)
#if defined(__GNUC__)
#define IRAFL_API __attribute__ ((dllexport))
#else
#define IRAFL_API __declspec(dllexport)
#endif
#else
#if defined(__GNUC__)
#define IRAFL_API __attribute__ ((dllimport))
#else
#define IRAFL_API __declspec(dllimport)
#endif
#endif
#define IRAFL_LOCAL
#else
#if __GNUC__ >= 4
#define IRAFL_API __attribute__ ((visibility ("default")))
#define IRAFL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define IRAFL_API
#define IRAFL_LOCAL
#endif
#endif

#if !defined(IRAFL_STRING_CONSTANCE)
#define IRAFL_STRING_CONSTANCE 1
#define IRAFL_VERSION_MAJORITY 1
#define IRAFL_VERSION_MINORITY 0
#define IRAFL_VERSION_APPENDED 0
#define IRAFL_VERSION_STR(s) #s
#define IRAFL_VERSION_DEC(majority, minority, appended) IRAFL_VERSION_STR(majority) "." IRAFL_VERSION_STR(minority) "." IRAFL_VERSION_STR(appended)
#define IRAFL_VERSION IRAFL_VERSION_DEC(IRAFL_VERSION_MAJORITY, IRAFL_VERSION_MINORITY, IRAFL_VERSION_APPENDED)
#define IRAFL_INFO \
"(I)mage (R)endering (A)nd (F)ilter (L)ibrary\n" \
"Version: " IRAFL_VERSION
#endif

#if !defined(IRAFL_FILTER_TYPE)
#define IRAFL_FILTER_TYPE 1
#define IRAFL_FILTER_TYPE_GRAY      1
#define IRAFL_FILTER_TYPE_GAUSSBLUR 2
#define IRAFL_FILTER_TYPE_LUT3D     4
#define IRAFL_FILTER_TYPE_VORONOI   8
#endif

typedef struct _VP
{
    double r;
    uint32_t k;
} VP;

#endif
