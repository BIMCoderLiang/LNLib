/*
 * Author:
 * 2023/10/17 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "KnotVectorUtils.h"
#include "UV.h"
#include "MathUtils.h"
#include "Polynomials.h"
#include "ValidationUtils.h"
#include "LNLibExceptions.h"
#include <algorithm>

namespace LNLib
{
	void InsertMidKnotCore(std::vector<double>& unqiueKnotVector, std::vector<double>& insert, int limitNumber)
	{
		if (insert.size() == limitNumber)
		{
			return;
		}
		else
		{
			double standard = Constants::DoubleEpsilon;
			int index = -1;
			for (int i = 0; i < unqiueKnotVector.size() - 1; i++)
			{
				double delta = unqiueKnotVector[i + 1] - unqiueKnotVector[i];
				if (MathUtils::IsGreaterThan(delta, standard))
				{
					standard = delta;
					index = i;
				}
			}
			double current = unqiueKnotVector[index] + standard / 2.0;
			unqiueKnotVector.emplace_back(current);
			std::sort(unqiueKnotVector.begin(), unqiueKnotVector.end());
			insert.emplace_back(current);

			InsertMidKnotCore(unqiueKnotVector, insert, limitNumber);
		}
	}
}

int LNLib::KnotVectorUtils::GetContinuity(int degree, const std::vector<double>& knotVector, double knot)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must be greater than zero.");
	int multi = Polynomials::GetKnotMultiplicity(knotVector, knot);
	return degree - multi;
}

std::vector<double> LNLib::KnotVectorUtils::Rescale(const std::vector<double>& knotVector, double min, double max)
{
	double originMin = knotVector[0];
	double originMax = knotVector[knotVector.size() - 1];
	double k = (max - min) / (originMax - originMin);

	int size = knotVector.size();
	std::vector<double> result(size);
	for (int i = 0; i < size; i++)
	{
		result[i] = k * (knotVector[i] - originMin) + min;
	}
	return result;
}

std::vector<double> LNLib::KnotVectorUtils::GetInsertedKnotElement(
	int degree,
	const std::vector<double>& knotVector,
	double startParam,
	double endParam)
{
	VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must be >= 0.");
	VALIDATE_ARGUMENT(!knotVector.empty(), "knotVector", "KnotVector must not be empty.");
	VALIDATE_ARGUMENT(startParam <= endParam, "startParam/endParam", "startParam must <= endParam.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "Must be non-decreasing.");

	std::vector<double> result;

	int startMulti = 0;
	for (double k : knotVector) {
		if (MathUtils::IsAlmostEqualTo(k, startParam)) {
			++startMulti;
		}
		else {
			break; 
		}
	}

	if (startMulti < degree + 1) { 
		int need = (degree + 1) - startMulti;
		for (int i = 0; i < need; ++i) {
			result.push_back(startParam);
		}
	}

	int endMulti = 0;
	for (auto it = knotVector.rbegin(); it != knotVector.rend(); ++it) {
		if (MathUtils::IsAlmostEqualTo(*it, endParam)) {
			++endMulti;
		}
		else {
			break;
		}
	}

	if (endMulti < degree + 1) {
		int need = (degree + 1) - endMulti;
		for (int i = 0; i < need; ++i) {
			result.push_back(endParam);
		}
	}

	return result;
}

std::map<double, int> LNLib::KnotVectorUtils::GetKnotMultiplicityMap(const std::vector<double>& knotVector)
{
	std::map<double, int> result;

	for (double knot : knotVector) {

		++result[knot];
	}

	return result;
}

std::map<double, int> LNLib::KnotVectorUtils::GetInternalKnotMultiplicityMap(const std::vector<double>& knotVector)
{
	auto result = GetKnotMultiplicityMap(knotVector);
	if (!result.empty()) 
	{
		result.erase(result.begin());
		auto lastElementIterator = prev(result.end(), 1);
		result.erase(lastElementIterator);
	}
	return result;
}

void LNLib::KnotVectorUtils::GetInsertedKnotElement(const std::vector<double>& knotVector0, const std::vector<double>& knotVector1, std::vector<double>& insertElements0, std::vector<double>& insertElements1)
{
	std::map<double, int> map0 = GetKnotMultiplicityMap(knotVector0);
	std::map<double, int> map1 = GetKnotMultiplicityMap(knotVector1);

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

std::vector<std::vector<double>> LNLib::KnotVectorUtils::GetInsertedKnotElements(const std::vector<std::vector<double>>& knotVectors)
{
	std::map<double, int> finalMap;
	for (int i = 0; i < knotVectors.size(); i++)
	{
		auto kv = knotVectors[i];
		std::map<double, int> map = GetKnotMultiplicityMap(kv);
		for (auto it = map.begin(); it != map.end(); ++it)
		{
			double key = it->first;
			int count = it->second;
			auto got = finalMap.find(key);
			if (got == finalMap.end())
			{
				finalMap.insert({ key, count });
			}
			else
			{
				int currentCount = got->second;
				if (currentCount < count)
				{
					finalMap[key] = count;
				}
			}
		}
	}

	std::vector<std::vector<double>> result;
	for (int i = 0; i < knotVectors.size(); i++)
	{
		auto kv = knotVectors[i];
		std::map<double, int> map = GetKnotMultiplicityMap(kv);

		std::vector<double> insertElements;
		for (auto it = finalMap.begin(); it != finalMap.end(); ++it)
		{
			double key = it->first;
			int count = it->second;

			auto got = map.find(key);
			if (got == map.end())
			{
				for (int j = 0; j < count; j++)
				{
					insertElements.emplace_back(key);
				}
			}
			else
			{
				int currentCount = got->second;
				int times = count - currentCount;
				for (int j = 0; j < times; j++)
				{
					insertElements.emplace_back(key);
				}
			}
		}
		result.emplace_back(insertElements);
	}
	return result;
}

std::vector<double> LNLib::KnotVectorUtils::GetInsertedKnotElements(int insertKnotsNumber, const std::vector<double>& knotVector)
{
	std::vector<double> unqiueKnotVector = knotVector;
	unqiueKnotVector.erase(std::unique(unqiueKnotVector.begin(), unqiueKnotVector.end()), unqiueKnotVector.end());

	std::vector<double> insert;
	InsertMidKnotCore(unqiueKnotVector, insert, insertKnotsNumber);
	std::sort(insert.begin(), insert.end());
	return insert;
}

bool LNLib::KnotVectorUtils::IsUniform(const std::vector<double>& knotVector)
{
	auto map = GetKnotMultiplicityMap(knotVector);
	if (map.empty())
	{
		return false;
	}
	std::vector<double> knots;
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		double key = it->first;
		knots.emplace_back(key);
	}
	std::sort(knots.begin(), knots.end());
	if (map[knots[0]] != map[knots[knots.size() - 1]])
	{
		return false;
	}
	double standard = knots[1] - knots[0];
	for (int i = 1; i < knots.size() - 1; i++)
	{
		double current = knots[i];
		double next = knots[i + 1];

		double gap = next - current;
		if (!MathUtils::IsAlmostEqualTo(gap, standard))
		{
			return false;
		}
		int currentCounts = map[current];
		int nextCounts = map[next];
		if (currentCounts != nextCounts)
		{
			if (i + 1 == knots.size() - 1)
			{
				continue;
			}
			else
			{
				return false;
			}
		}
	}
	return true;
}



