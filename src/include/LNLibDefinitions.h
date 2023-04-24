#pragma once

#define DLL_EXPORT __declspec(dllexport)
#define DLL_IMPORT __declspec(dllimport)

#ifdef LNLIB_HOME
#define LNLIB_EXPORT DLL_EXPORT
#else
#define LNLIB_EXPORT DLL_IMPORT
#endif


