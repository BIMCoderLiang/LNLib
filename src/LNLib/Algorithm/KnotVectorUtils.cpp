/*
 * Author:
 * 2023/10/17 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "KnotVectorUtils.h"
#include "UV.h"
#include "MathUtils.h"
#include "Polynomials.h"
#include "ValidationUtils.h"
#include "LNLibExceptions.h"
#include <algorithm>

using namespace LNLib;

int LNLib::KnotVectorUtils::GetContinuity(int degree, const std::vector<double>& knotVector, double knot)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	int multi = Polynomials::GetKnotMultiplicity(knotVector, knot);
	return degree - multi;
}

std::vector<double> LNLib::KnotVectorUtils::GetInsertedKnotElement(int degree, const std::vector<double>& knotVector, double startParam, double endParam)
{
	VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must greater than or equals zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");

	std::vector<double> result;
	int startMulti = Polynomials::GetKnotMultiplicity(knotVector, startParam);
	if (startMulti < degree)
	{
		for (int i = 0; i < degree - startMulti; i++)
		{
			result.emplace_back(startParam);
		}
	}

	int endMulti = Polynomials::GetKnotMultiplicity(knotVector, endParam);
	if (endMulti < degree)
	{
		for (int i = 0; i < degree - endMulti; i++)
		{
			result.emplace_back(endParam);
		}
	}

	return result;
}

std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> LNLib::KnotVectorUtils::GetKnotMultiplicityMap(const std::vector<double>& knotVector)
{
	std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> result;

	for (int i = 0; i < knotVector.size(); i++)
	{
		auto got = result.find(knotVector[i]);
		if (got == result.end())
		{
			int multi = Polynomials::GetKnotMultiplicity(knotVector, knotVector[i]);
			result.insert({ knotVector[i], multi });
		}
	}
	return result;
}

void LNLib::KnotVectorUtils::GetInsertedKnotElement(const std::vector<double> knotVector0, const std::vector<double> knotVector1, std::vector<double>& insertElements0, std::vector<double>& insertElements1)
{
	std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> map0 = GetKnotMultiplicityMap(knotVector0);
	std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> map1 = GetKnotMultiplicityMap(knotVector1);

	for (auto it = map0.begin(); it != map0.end(); ++it)
	{
		double key0 = it->first;
		int count0 = it->second;

		auto got = map1.find(key0);
		if (got == map1.end())
		{
			for (int i = 0; i < count0; i++)
			{
				insertElements1.emplace_back(key0);
			}
		}
		else
		{
			int count1 = got->second;
			if (count0 > count1)
			{
				int times = count0 - count1;
				for (int j = 0; j < times; j++)
				{
					insertElements1.emplace_back(key0);
				}
			}
			else
			{
				int times = count1 - count0;
				for (int j = 0; j < times; j++)
				{
					insertElements0.emplace_back(key0);
				}
			}
		}
	}

	for (auto it = map1.begin(); it != map1.end(); ++it)
	{
		double key1 = it->first;
		int count1 = it->second;

		auto got = map0.find(key1);
		if (got == map0.end())
		{
			for (int i = 0; i < count1; i++)
			{
				insertElements0.emplace_back(key1);
			}
		}
	}

	std::sort(insertElements0.begin(), insertElements0.end());
	std::sort(insertElements1.begin(), insertElements1.end());
}


