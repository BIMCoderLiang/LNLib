/*
 * Author:
 * 2024/01/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Integrator.h"
#include "FFT.h"
#include "Constants.h"

#include <cmath>

namespace LNLib
{
    double Integrator::Simpson(IntegrationFunction& function, void* customData, double start, double end)
    {
        double st = (function)(start, customData);
        double mt = (function)((start + end) / 2.0, customData);
        double et = (function)((end), customData);
        double result = ((end - start) / 6.0) * (st + 4 * mt + et);
        return result;
    }

	double Integrator::Simpson(BinaryIntegrationFunction& function, void* customData, 
				double u1, double u2, double v1, double v2)
	{
        double du = u2 - u1;
        double dv = v2 - v1;
        double hdu = 0.5 * du;
        double hdv = 0.5 * dv;
		
        // Sample 9 points with weights
        double uvw[] = {
            u1,       v1,       1,
            u1,       v1 + hdv, 4,
            u1,       v2,       1,
            u1 + hdu, v1,       4,
            u1 + hdu, v1 + hdv, 16,
            u1 + hdu, v2,       4,
            u2,       v1,       1,
            u2,       v1 + hdv, 4,
            u2,       v2,       1,
            };
        
        double sum = 0;
        for(int i=0;i<9;++i)
        {
            double* base = uvw + i*3;
            double u = base[0];
            double v = base[1];
            double w = base[2];
            double f = function(u, v, customData);
            sum += w * f;
        }

        sum *= du * dv / 36;
        return sum;
	}

    const std::vector<double> Integrator::GaussLegendreAbscissae =
    {
        -0.0640568928626056260850430826247450385909,
        0.0640568928626056260850430826247450385909,
        -0.1911188674736163091586398207570696318404,
        0.1911188674736163091586398207570696318404,
        -0.3150426796961633743867932913198102407864,
        0.3150426796961633743867932913198102407864,
        -0.4337935076260451384870842319133497124524,
        0.4337935076260451384870842319133497124524,
        -0.5454214713888395356583756172183723700107,
        0.5454214713888395356583756172183723700107,
        -0.6480936519369755692524957869107476266696,
        0.6480936519369755692524957869107476266696,
        -0.7401241915785543642438281030999784255232,
        0.7401241915785543642438281030999784255232,
        -0.8200019859739029219539498726697452080761,
        0.8200019859739029219539498726697452080761,
        -0.8864155270044010342131543419821967550873,
        0.8864155270044010342131543419821967550873,
        -0.9382745520027327585236490017087214496548,
        0.9382745520027327585236490017087214496548,
        -0.9747285559713094981983919930081690617411,
        0.9747285559713094981983919930081690617411,
        -0.9951872199970213601799974097007368118745,
        0.9951872199970213601799974097007368118745,
    };

    const std::vector<double> Integrator::GaussLegendreWeights =
    {
        0.1279381953467521569740561652246953718517,
        0.1279381953467521569740561652246953718517,
        0.1258374563468282961213753825111836887264,
        0.1258374563468282961213753825111836887264,
        0.121670472927803391204463153476262425607,
        0.121670472927803391204463153476262425607,
        0.1155056680537256013533444839067835598622,
        0.1155056680537256013533444839067835598622,
        0.1074442701159656347825773424466062227946,
        0.1074442701159656347825773424466062227946,
        0.0976186521041138882698806644642471544279,
        0.0976186521041138882698806644642471544279,
        0.086190161531953275917185202983742667185,
        0.086190161531953275917185202983742667185,
        0.0733464814110803057340336152531165181193,
        0.0733464814110803057340336152531165181193,
        0.0592985849154367807463677585001085845412,
        0.0592985849154367807463677585001085845412,
        0.0442774388174198061686027482113382288593,
        0.0442774388174198061686027482113382288593,
        0.0285313886289336631813078159518782864491,
        0.0285313886289336631813078159518782864491,
        0.0123412297999871995468056670700372915759,
        0.0123412297999871995468056670700372915759,
    };

    std::vector<double> Integrator::ChebyshevSeries(int size)
    {
        std::vector<double> series(size);

        int lenw = series.size() - 1;
        int j, k, l, m;
        double cos2, sin1, sin2, hl;

        cos2 = 0;
        sin1 = 1;
        sin2 = 1;
        hl = 0.5;
        k = lenw;
        l = 2;
        while (l < k - l - 1) 
        {
            series[0] = hl * 0.5;
            for (j = 1; j <= l; j++) 
            {
                series[j] = hl / (1 - 4 * j * j);
            }
            series[l] *= 0.5;
            dfct(l, 0.5 * cos2, sin1, series);
            cos2 = std::sqrt(2 + cos2);
            sin1 /= cos2;
            sin2 /= 2 + cos2;
            series[k] = sin2;
            series[k - 1] = series[0];
            series[k - 2] = series[l];
            k -= 3;
            m = l;
            while (m > 1) 
            {
                m >>= 1;
                for (j = m; j <= l - m; j += (m << 1)) 
                {
                    series[k] = series[j];
                    k--;
                }
            }
            hl *= 0.5;
            l *= 2;
        }
        return series;
    }
    double Integrator::ClenshawCurtisQuadrature(IntegrationFunction& function, void* customData, double start, double end, std::vector<double>& series, double epsilon)
    {
        double integration;
        int j, k, l;
        double err, esf, eref, erefh, hh, ir, iback, irback, ba, ss, x, y, fx, errir;
        int lenw = series.size() - 1;
        esf = 10;
        ba = 0.5 * (end - start);
        ss = 2 * series[lenw];
        x = ba * series[lenw];
        series[0] = 0.5 * (function)(start, customData);
        series[3] = 0.5 * (function)(end, customData);
        series[2] = (function)(start + x, customData);
        series[4] = (function)(end - x, customData);
        series[1] = (function)(start + ba, customData);
        eref = 0.5 * (fabs(series[0]) + std::fabs(series[1]) + std::fabs(series[2]) + std::fabs(series[3]) + std::fabs(series[4]));
        series[0] += series[3];
        series[2] += series[4];
        ir = series[0] + series[1] + series[2];
        integration = series[0] * series[lenw - 1] + series[1] * series[lenw - 2] + series[2] * series[lenw - 3];
        erefh = eref * std::sqrt(epsilon);
        eref *= epsilon;
        hh = 0.25;
        l = 2;
        k = lenw - 5;
        do {
            iback = integration;
            irback = ir;
            x = ba * series[k + 1];
            y = 0;
            integration = series[0] * series[k];
            for (j = 1; j <= l; j++) {
                x += y;
                y += ss * (ba - x);
                fx = (function)(start + x, customData) + (function)(end - x, customData);
                ir += fx;
                integration += series[j] * series[k - j] + fx * series[k - j - l];
                series[j + l] = fx;
            }
            ss = 2 * series[k + 1];
            err = esf * l * std::fabs(integration - iback);
            hh *= 0.25;
            errir = hh * std::fabs(ir - 2 * irback);
            l *= 2;
            k -= l + 2;
        } while ((err > erefh || errir > eref) && k > 4 * l);
        integration *= end - start;
        if (err > erefh || errir > eref)
        {
            err *= -fabs(end - start);
        }
        else
        {
            err = eref * std::fabs(end - start);
        }
        return integration;
    }
    double Integrator::ClenshawCurtisQuadrature2(IntegrationFunction& function, void* customData, double start, double end, std::vector<double> series, double epsilon)
    {
        double integration;
        int j, k, l;
        double err, esf, eref, erefh, hh, ir, iback, irback, ba, ss, x, y, fx, errir;
        int lenw = series.size() - 1;
        esf = 10;
        ba = 0.5 * (end - start);
        ss = 2 * series[lenw];
        x = ba * series[lenw];
        series[0] = 0.5 * (function)(start, customData);
        series[3] = 0.5 * (function)(end, customData);
        series[2] = (function)(start + x, customData);
        series[4] = (function)(end - x, customData);
        series[1] = (function)(start + ba, customData);
        eref = 0.5 * (fabs(series[0]) + std::fabs(series[1]) + std::fabs(series[2]) + std::fabs(series[3]) + std::fabs(series[4]));
        series[0] += series[3];
        series[2] += series[4];
        ir = series[0] + series[1] + series[2];
        integration = series[0] * series[lenw - 1] + series[1] * series[lenw - 2] + series[2] * series[lenw - 3];
        erefh = eref * std::sqrt(epsilon);
        eref *= epsilon;
        hh = 0.25;
        l = 2;
        k = lenw - 5;
        do {
            iback = integration;
            irback = ir;
            x = ba * series[k + 1];
            y = 0;
            integration = series[0] * series[k];
            for (j = 1; j <= l; j++) {
                x += y;
                y += ss * (ba - x);
                fx = (function)(start + x, customData) + (function)(end - x, customData);
                ir += fx;
                integration += series[j] * series[k - j] + fx * series[k - j - l];
                series[j + l] = fx;
            }
            ss = 2 * series[k + 1];
            err = esf * l * std::fabs(integration - iback);
            hh *= 0.25;
            errir = hh * std::fabs(ir - 2 * irback);
            l *= 2;
            k -= l + 2;
        } while ((err > erefh || errir > eref) && k > 4 * l);
        integration *= end - start;
        if (err > erefh || errir > eref)
        {
            err *= -fabs(end - start);
        }
        else
        {
            err = eref * std::fabs(end - start);
        }
        return integration;
    }
}


