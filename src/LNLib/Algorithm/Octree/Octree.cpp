/*
 * Author:
 * 2026/06/26 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Octree.h"
#include "OctreeCollector.h"
#include "XYZ.h"
#include "MathUtils.h"

#include <cmath>
#include <algorithm>
#include <limits>


LNLib::Octant::Octant(LNLib::XYZ center, double extent, std::vector<int> indices, bool isLeaf)
    : Center(center), Extent(extent), PointIndices(std::move(indices)), IsLeaf(isLeaf) {
    Children.resize(8, nullptr);
}

LNLib::Octant::~Octant() {
    for (auto* child : Children) {
        delete child;
    }
}

bool LNLib::Inside(const LNLib::XYZ& query, double radius, LNLib::Octant* octant) {
    LNLib::XYZ queryOffset = query.Substract(octant->Center);
    double ax = std::fabs(queryOffset.X());
    double ay = std::fabs(queryOffset.Y());
    double az = std::fabs(queryOffset.Z());
    return (MathUtils::IsLessThan(ax + radius,octant->Extent)) && (MathUtils::IsLessThan(ay + radius,octant->Extent)) && (MathUtils::IsLessThan(az + radius,octant->Extent));
}

bool LNLib::Overlaps(const LNLib::XYZ& query, double radius, LNLib::Octant* octant) {
    LNLib::XYZ queryOffset = query.Substract(octant->Center);
    double ax = std::fabs(queryOffset.X());
    double ay = std::fabs(queryOffset.Y());
    double az = std::fabs(queryOffset.Z());

    double maxDist = radius + octant->Extent;
    if (MathUtils::IsGreaterThan(ax,maxDist) || MathUtils::IsGreaterThan(ay,maxDist) || MathUtils::IsGreaterThan(az,maxDist)) return false;

    int count = 0;
    if (MathUtils::IsLessThan(ax,octant->Extent)) count++;
    if (MathUtils::IsLessThan(ay,octant->Extent)) count++;
    if (MathUtils::IsLessThan(az,octant->Extent)) count++;
    if (count >= 2) return true;

    double xDiff = std::max(ax - octant->Extent, 0.0);
    double yDiff = std::max(ay - octant->Extent, 0.0);
    double zDiff = std::max(az - octant->Extent, 0.0);

    LNLib::XYZ diff(xDiff, yDiff, zDiff);
    return  MathUtils::IsLessThan(diff.SqrLength(),(radius * radius));
}

bool LNLib::Contains(const LNLib::XYZ& query, double radius, LNLib::Octant* octant) {
    LNLib::XYZ queryOffset = query.Substract(octant->Center);
    double ax = std::fabs(queryOffset.X());
    double ay = std::fabs(queryOffset.Y());
    double az = std::fabs(queryOffset.Z());
    LNLib::XYZ farthestCorner(ax + octant->Extent, ay + octant->Extent, az + octant->Extent);
    return farthestCorner.Length() < radius;
}

LNLib::Octant* OctreeRecursiveBuild(LNLib::Octant* root, const std::vector<LNLib::XYZ>& db,
    const LNLib::XYZ& center, double extent,
    const std::vector<int>& pointIndices,
    int leafSize, double minExtent) {
    if (pointIndices.empty()) return nullptr;

    if (root == nullptr) {
        root = new LNLib::Octant(center, extent, pointIndices, true);
    }

    if ((int)pointIndices.size() <= leafSize || extent <= minExtent) {
        root->IsLeaf = true;
    }
    else {
        root->IsLeaf = false;
        std::vector<std::vector<int>> childrenPointIndices(8);

        for (int pointIdx : pointIndices) {
            const LNLib::XYZ& pointDb = db[pointIdx];
            int mortonCode = 0;
            if (LNLib::MathUtils::IsGreaterThan(pointDb.X(),center.X())) mortonCode |= 1;
            if (LNLib::MathUtils::IsGreaterThan(pointDb.Y(),center.Y())) mortonCode |= 2;
            if (LNLib::MathUtils::IsGreaterThan(pointDb.Z(),center.Z())) mortonCode |= 4;
            childrenPointIndices[mortonCode].push_back(pointIdx);
        }

        double factor[2] = { -0.5, 0.5 };
        for (int i = 0; i < 8; ++i) {
            LNLib::XYZ childCenter(
                center.X() + factor[(i & 1) > 0] * extent,
                center.Y() + factor[(i & 2) > 0] * extent,
                center.Z() + factor[(i & 4) > 0] * extent
            );
            double childExtent = 0.5 * extent;
            root->Children[i] = OctreeRecursiveBuild(root->Children[i], db,
                childCenter, childExtent,
                childrenPointIndices[i],
                leafSize, minExtent);
        }
    }
    return root;
}

LNLib::Octant* LNLib::OctreeConstruction(const std::vector<LNLib::XYZ>& db, int leafSize, double minExtent) {
    int n = (int)db.size();
    if (n == 0) return nullptr;

    LNLib::XYZ dbMin(std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max());
    LNLib::XYZ dbMax(std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest());

    for (const auto& p : db) {

        if (LNLib::MathUtils::IsLessThan(p.X(),dbMin.X())) dbMin.X() = p.X();
        if (LNLib::MathUtils::IsLessThan(p.Y(),dbMin.Y())) dbMin.Y() = p.Y();
        if (LNLib::MathUtils::IsLessThan(p.Z(),dbMin.Z())) dbMin.Z() = p.Z();

        if (LNLib::MathUtils::IsGreaterThan(p.X(),dbMax.X())) dbMax.X() = p.X();
        if (LNLib::MathUtils::IsGreaterThan(p.Y(),dbMax.Y())) dbMax.Y() = p.Y();
        if (LNLib::MathUtils::IsGreaterThan(p.Z(),dbMax.Z())) dbMax.Z() = p.Z();
    }

    double dbExtent = std::max({ dbMax.X() - dbMin.X(),
                                 dbMax.Y() - dbMin.Y(),
                                 dbMax.Z() - dbMin.Z() }) * 0.5;

    LNLib::XYZ dbCenter(dbMin.X() + dbExtent,
        dbMin.Y() + dbExtent,
        dbMin.Z() + dbExtent);

    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;

    return OctreeRecursiveBuild(nullptr, db, dbCenter, dbExtent, indices, leafSize, minExtent);
}

bool LNLib::OctreeRadiusSearchFast(LNLib::Octant* root, const std::vector<LNLib::XYZ>& db,
    RadiusNnResultSet& resultSet, const LNLib::XYZ& query) {
    if (root == nullptr) return false;

    if (Contains(query, resultSet.WorstDist(), root)) {
        for (int idx : root->PointIndices) {
            double dist = query.Distance(db[idx]);
            resultSet.AddPoint(dist, idx);
        }
        return false;
    }

    if (root->IsLeaf && !root->PointIndices.empty()) {
        for (int idx : root->PointIndices) {
            double dist = query.Distance(db[idx]);
            resultSet.AddPoint(dist, idx);
        }
        return Inside(query, resultSet.WorstDist(), root);
    }

    for (LNLib::Octant* child : root->Children) {
        if (child == nullptr) continue;
        if (!Overlaps(query, resultSet.WorstDist(), child)) continue;
        if (LNLib::OctreeRadiusSearchFast(child, db, resultSet, query)) return true;
    }
    return LNLib::Inside(query, resultSet.WorstDist(), root);
}

bool LNLib::OctreeRadiusSearch(LNLib::Octant* root, const std::vector<LNLib::XYZ>& db,
    RadiusNnResultSet& resultSet, const LNLib::XYZ& query) {
    if (root == nullptr) return false;

    if (root->IsLeaf && !root->PointIndices.empty()) {
        for (int idx : root->PointIndices) {
            double dist = query.Distance(db[idx]);
            resultSet.AddPoint(dist, idx);
        }
        return Inside(query, resultSet.WorstDist(), root);
    }

    int mortonCode = 0;
    if (LNLib::MathUtils::IsGreaterThan(query.X(),root->Center.X())) mortonCode |= 1;
    if (LNLib::MathUtils::IsGreaterThan(query.Y(),root->Center.Y())) mortonCode |= 2;
    if (LNLib::MathUtils::IsGreaterThan(query.Z(),root->Center.Z())) mortonCode |= 4;

    if (LNLib::OctreeRadiusSearch(root->Children[mortonCode], db, resultSet, query))
        return true;

    for (int c = 0; c < 8; ++c) {
        if (c == mortonCode || root->Children[c] == nullptr) continue;
        if (!Overlaps(query, resultSet.WorstDist(), root->Children[c])) continue;
        if (LNLib::OctreeRadiusSearch(root->Children[c], db, resultSet, query))
            return true;
    }
    return Inside(query, resultSet.WorstDist(), root);
}

bool LNLib::OctreeKnnSearch(LNLib::Octant* root, const std::vector<LNLib::XYZ>& db,
    KnnResultSet& resultSet, const LNLib::XYZ& query) {
    if (root == nullptr) return false;

    if (root->IsLeaf && !root->PointIndices.empty()) {
        for (int idx : root->PointIndices) {
            double dist = query.Distance(db[idx]);
            resultSet.AddPoint(dist, idx);
        }
        return Inside(query, resultSet.WorstDist(), root);
    }

    int mortonCode = 0;
    if (LNLib::MathUtils::IsGreaterThan(query.X(),root->Center.X())) mortonCode |= 1;
    if (LNLib::MathUtils::IsGreaterThan(query.Y(),root->Center.Y())) mortonCode |= 2;
    if (LNLib::MathUtils::IsGreaterThan(query.Z(),root->Center.Z())) mortonCode |= 4;

    if (LNLib::OctreeKnnSearch(root->Children[mortonCode], db, resultSet, query))
        return true;

    for (int c = 0; c < 8; ++c) {
        if (c == mortonCode || root->Children[c] == nullptr) continue;
        if (!Overlaps(query, resultSet.WorstDist(), root->Children[c])) continue;
        if (LNLib::OctreeKnnSearch(root->Children[c], db, resultSet, query))
            return true;
    }
    return Inside(query, resultSet.WorstDist(), root);
}
