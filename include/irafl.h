#if !defined(_IRAFL_H_)
#define _IRAFL_H_

#include <cstdint>
#include <string>
#include "irafldef.h"

class IRAFL_API IRAFL
{
protected:
	class InnerImpl;
protected:
	IRAFL();
public:
	static const char* version;
	static const char* info;
public:
	virtual ~IRAFL();
	static IRAFL* Handle();
	virtual void Initialize();
	virtual void Terminate();
public:
	static const char* Version();
	static const char* Info();

	static void Execute(unsigned char*& dst, const unsigned char* src, const long& w, const long& h, const uint32_t& id, const void* params = nullptr);
	static const std::string& GetFilterName(const uint32_t& id);
protected:
	static IRAFL* instance;
protected:
	InnerImpl* inner = nullptr;
};
#define IRAFLH IRAFL::Handle()

#endif