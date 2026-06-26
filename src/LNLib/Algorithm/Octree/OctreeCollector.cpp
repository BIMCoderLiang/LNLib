/*
 * Author:
 * 2026/06/26 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "OctreeCollector.h"
#include "MathUtils.h"

LNLib::KnnResultSet::KnnResultSet(int capacity)
    : Capacity(capacity), Count(0), WorstDistValue(1e10), ComparisonCounter(0) {
    DistIndexList.resize(Capacity);
    for (int i = 0; i < Capacity; ++i) {
        DistIndexList[i] = DistIndex(WorstDistValue, 0);
    }
}

int LNLib::KnnResultSet::Size() const { return Count; }
bool LNLib::KnnResultSet::Full() const { return Count == Capacity; }
double LNLib::KnnResultSet::WorstDist() const { return WorstDistValue; }

void LNLib::KnnResultSet::AddPoint(double dist, int index) {
    ComparisonCounter++;

    if (MathUtils::IsGreaterThan(dist,WorstDistValue)) return;

    int i;
    if (Count < Capacity) {
        Count++;
        i = Count - 1;
    }
    else {
        i = Capacity - 1;
    }

    while (i > 0) {
        if (MathUtils::IsGreaterThan(DistIndexList[i - 1].Distance,dist)) {
            DistIndexList[i] = DistIndexList[i - 1];
            i--;
        }
        else {
            break;
        }
    }

    DistIndexList[i].Distance = dist;
    DistIndexList[i].Index = index;
    WorstDistValue = DistIndexList[Capacity - 1].Distance;
}

LNLib::RadiusNnResultSet::RadiusNnResultSet(double radius)
    : Radius(radius), Count(0), ComparisonCounter(0) {
}

int LNLib::RadiusNnResultSet::Size() const { return Count; }
double LNLib::RadiusNnResultSet::WorstDist() const { return Radius; }

void LNLib::RadiusNnResultSet::AddPoint(double dist, int index) {
    ComparisonCounter++;
    if (MathUtils::IsGreaterThan(dist,Radius)) return;

    Count++;
    DistIndexList.push_back(DistIndex(dist, index));
}




