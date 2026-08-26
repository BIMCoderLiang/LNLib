/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once

#define DLL_EXPORT __declspec(dllexport)
#define DLL_IMPORT __declspec(dllimport)

#if defined(WIN64) || defined(_WIN64) || defined(__WIN64__) || defined(__CYGWIN__)
    #if defined(LNLIB_HOME) && !defined(LNLIB_STATIC)
        #define LNLIB_EXPORT __declspec(dllexport)
    #elif !defined(LNLIB_HOME)
        #define LNLIB_EXPORT __declspec(dllimport)
    #else
        #define LNLIB_EXPORT
    #endif
#else
    /*
     * POSIX: no export attribute needed. CMake does not enable
     * -fvisibility=hidden, so all symbols are visible by default, and GCC
     * warns (-Wattributes) when visibility is applied to a typedef (used for
     * LN_NurbsCurve etc. in LNObject.h).
     */
    #define LNLIB_EXPORT
#endif
