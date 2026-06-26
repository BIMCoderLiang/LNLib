/*
 * Author:
 * 2026/06/26 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "XYZ.h"
#include "OctreeCollector.h"
#include <vector>

namespace LNLib
{
    class Octant {
    public:
        std::vector<Octant*> Children;
        LNLib::XYZ Center;
        double Extent;
        std::vector<int> PointIndices;
        bool IsLeaf;

        Octant(LNLib::XYZ center, double extent, std::vector<int> indices, bool isLeaf);
        ~Octant();
    };


    bool Inside(const LNLib::XYZ& query, double radius, Octant* octant);
    bool Overlaps(const LNLib::XYZ& query, double radius, Octant* octant);
    bool Contains(const LNLib::XYZ& query, double radius, Octant* octant);


    Octant* OctreeConstruction(const std::vector<LNLib::XYZ>& db, int leafSize, double minExtent);

    bool OctreeRadiusSearchFast(Octant* root, const std::vector<LNLib::XYZ>& db,
        RadiusNnResultSet& resultSet, const LNLib::XYZ& query);

    bool OctreeRadiusSearch(Octant* root, const std::vector<LNLib::XYZ>& db,
        RadiusNnResultSet& resultSet, const LNLib::XYZ& query);

    bool OctreeKnnSearch(Octant* root, const std::vector<LNLib::XYZ>& db,
        KnnResultSet& resultSet, const LNLib::XYZ& query);
}