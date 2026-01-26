/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <numeric>
#include <vector>
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/Point.h"
#include "eckit/types/FloatCompare.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

using ConvexSphericalPolygon = util::ConvexSphericalPolygon;

const double relative_error           = 0.00005;
const double recipRoot2               = 1 / sqrt(2);
const double recipRoot3               = 1 / sqrt(3);


CASE("test_convex_spherical_polygon_triag") {
    std::array<PointXYZ, 3> testTriangleVertices{PointXYZ{0, 1, 0}, PointXYZ{recipRoot2, 0, recipRoot2},
                                                 PointXYZ{1, 0, 0}};

    util::ConvexSphericalPolygon testTriangle(testTriangleVertices.data(), testTriangleVertices.size());

    SECTION("test_edge_normals") {
        std::array<PointXYZ, 3> expectedEdgeNormals = {PointXYZ{recipRoot2, 0, -1 * recipRoot2},
                                                       PointXYZ{0, recipRoot2, 0}, PointXYZ{0, 0, 1}};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT(eckit::types::is_approximately_equal(expectedEdgeNormals[i][j],
                    testTriangle.edge_normals()[i][j]));
            }
        }
    }

    SECTION("test_weight_computation") {
        const size_t numberTestPoints = 10;

        std::array<PointXYZ, numberTestPoints> candidatePoints = {
            PointXYZ{0, 1, 0},                                            // Vertex 0
            PointXYZ{recipRoot2, 0, recipRoot2},                          // Vertex 1
            PointXYZ{1, 0, 0},                                            // Vertex 2
            PointXYZ{recipRoot3, recipRoot3, recipRoot3},                 // Inside
            PointXYZ{recipRoot2, recipRoot2, 0},                          // On edge
            PointXYZ{recipRoot2, recipRoot2, 0.000001},                   // Just inside edge
            PointXYZ{0.01, 0.9999, 0.01},                                 // Just inside vertex
            PointXYZ{0, recipRoot2, recipRoot2},                          // Outside
            PointXYZ{-1 * recipRoot3, -1 * recipRoot3, -1 * recipRoot3},  // Opposite Side
            PointXYZ{recipRoot2, recipRoot2, -0.000001}                   // Just outside edge
        };

        std::array<bool, numberTestPoints> isPointInside = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
        size_t expectedInside  = std::accumulate(std::begin(isPointInside), std::end(isPointInside), 0);
        size_t expectedOutside = numberTestPoints - expectedInside;

        std::array<std::vector<double>, numberTestPoints> candidateWeights = {
            std::vector<double>{1, 0, 0},            // Vertex 0
            std::vector<double>{0, 1, 0},            // Vertex 1
            std::vector<double>{0, 0, 1},            // Vertex 2
            std::vector<double>{0.5774, 0.8165, 0},  // Inside
            std::vector<double>{0.7071, 0, 0.7071},  // On edge
            std::vector<double>{0.7071, 0, 0.7071},  // Just inside edge
            std::vector<double>{0.9999, 0.0141, 0},  // Just inside vertex
            std::vector<double>{0, 0, 0},            // Outside
            std::vector<double>{0, 0, 0},            // Opposite Side
            std::vector<double>{0, 0, 0}             // Just outside edge
        };

        size_t pointsInside  = 0;
        size_t pointsOutside = 0;

        std::vector<double> polygonWeights(testTriangle.size());

        for (size_t i = 0; i < numberTestPoints; ++i) {
            if (testTriangle.compute_vertex_weights(candidatePoints[i], polygonWeights.data(), polygonWeights.size()) ==
                0) {
                EXPECT((isPointInside[i] == 0));
                pointsOutside += 1;
            }
            else {
                EXPECT((isPointInside[i] == 1));
                pointsInside += 1;
                for (size_t j = 0; j < 3; ++j) {
                    EXPECT(eckit::types::is_approximately_equal(polygonWeights[j], candidateWeights[i][j],
                                                                relative_error));
                }
            }
        }
        Log::info() << "Points in/out: " << pointsInside << "/" << pointsOutside << std::endl;
        EXPECT(pointsOutside == expectedOutside);
        EXPECT(pointsInside == expectedInside);
    }
}


CASE("test_spherical_polygon_nonplanar_quad") {
    std::array<PointXYZ, 4> testQuadVertices{PointXYZ{1,0,0},PointXYZ{recipRoot2,recipRoot2,0},{recipRoot3,recipRoot3,recipRoot3},{recipRoot2,0,recipRoot2}};

    util::ConvexSphericalPolygon testQuad(testQuadVertices.data(), testQuadVertices.size());
    const size_t numberTestPoints = 11;

    std::array<PointXYZ, numberTestPoints> candidatePoints = {
        PointXYZ{1, 0, 0},                             // Vertex 0
        PointXYZ{recipRoot2, recipRoot2, 0},           // Vertex 1
        PointXYZ{recipRoot3, recipRoot3, recipRoot3},  // Vertex2
        PointXYZ{recipRoot2, 0, recipRoot2},           // Vertex 3
        PointXYZ{0.8246, 0.4, 0.4},                    // Inside
        PointXYZ{0.5 * sqrt(3), 0.5, 0},            // On edge
        PointXYZ{0.5 * sqrt(3), 0.5, 0.000001},     // Just inside edge
        PointXYZ{0.5918, 0.57, 0.57},                  // Just inside vertex
        PointXYZ{recipRoot2, -1 * recipRoot2, 0},      // Outside
        PointXYZ{-0.8246, -0.4, -0.4},                 // Opposite side
        PointXYZ{0.5 * sqrt(3), 0.5, -0.000001},   // Just outside edge
    };

    std::array<bool, numberTestPoints> isPointInside = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
    const size_t expectedInside  = 8;
    const size_t expectedOutside                     = numberTestPoints - expectedInside;

    SECTION("test_edge_normals") {
        double recipRoot6 = recipRoot2 * recipRoot3;

        std::array<PointXYZ, 4> expectedEdgeNormals = {
            PointXYZ{0, 0, recipRoot2}, PointXYZ{recipRoot6, -1 * recipRoot6, 0},
            PointXYZ{recipRoot6, 0, -1 * recipRoot6}, PointXYZ{0, recipRoot2, 0}};
        
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 3; ++j) {
                EXPECT(
                    eckit::types::is_approximately_equal(expectedEdgeNormals[i][j], testQuad.edge_normals()[i][j]));
            }
        }
    }

    SECTION("test_weight_computation") {
        std::array<std::vector<double>, numberTestPoints> candidateWeights = {
            std::vector<double>{1, 0, 0, 0},
            std::vector<double>{0, 1, 0, 0},
            std::vector<double>{0, 0, 1, 0},
            std::vector<double>{0, 0, 0, 1},
            std::vector<double>{0.23847, 0.26325, 0.37043, 0.26325},
            std::vector<double>{0.3660, 0.7071, 0, 0},
            std::vector<double>{0.3660, 0.7071, 0, 0},
            std::vector<double>{0.0075, 0.0203, 0.9624, 0.0203},
            std::vector<double>{0, 0, 0, 0},
            std::vector<double>{0, 0, 0, 0},
            std::vector<double>{0, 0, 0, 0}  //
        };

        size_t pointsInside  = 0;
        size_t pointsOutside = 0;

        std::vector<double> polygonWeights(testQuad.size());

        for (size_t i = 0; i < numberTestPoints; ++i) {
            if (testQuad.compute_vertex_weights(candidatePoints[i], polygonWeights.data(), polygonWeights.size()) ==
                0) {
                EXPECT((isPointInside[i] == 0));
                pointsOutside += 1;
            }
            else {
                EXPECT((isPointInside[i] == 1));
                pointsInside += 1;
                for (size_t j = 0; j < 4; ++j) {
                    EXPECT(eckit::types::is_approximately_equal(polygonWeights[j], candidateWeights[i][j],
                                                                relative_error));
                }
            }
        }

        Log::info() << "Points in/out: " << pointsInside << "/" << pointsOutside << std::endl;
        EXPECT(pointsOutside == expectedOutside);
        EXPECT(pointsInside == expectedInside);

        // Test whether rotated polygon gives the same weights
        const std::array<PointXYZ, 4> testQuadRotatedVertices = {PointXYZ{recipRoot2, 0, recipRoot2}, PointXYZ{1, 0, 0},
                                                                 PointXYZ{recipRoot2, recipRoot2, 0},
                                                                 PointXYZ{recipRoot3, recipRoot3, recipRoot3}};

        util::ConvexSphericalPolygon testQuadRotated(testQuadRotatedVertices.data(), testQuadRotatedVertices.size());
        size_t pointsInsideRotated = 0;
        size_t pointsOutsideRotated = 0;

        std::vector<double> polygonWeightsRotated(testQuadRotated.size());

        for (size_t i = 0; i < numberTestPoints; ++i) {
            if (testQuadRotated.compute_vertex_weights(candidatePoints[i], polygonWeightsRotated.data(),
                                                       polygonWeightsRotated.size()) == 0) {
                EXPECT((isPointInside[i] == 0));
                pointsOutsideRotated += 1;
            }
            else {
                EXPECT((isPointInside[i] == 1));
                pointsInsideRotated += 1;
                for (size_t j = 0; j < 4; ++j) {
                    EXPECT(eckit::types::is_approximately_equal(
                        polygonWeightsRotated[j], candidateWeights[i][testQuadRotated.previous(j)], relative_error));
                }
            }
        }

        Log::info() << "Points in/out: " << pointsInsideRotated << "/" << pointsOutsideRotated << std::endl;
        EXPECT(pointsOutsideRotated == expectedOutside);
        EXPECT(pointsInsideRotated == expectedInside);
    }
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}