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
#include "atlas/interpolation/Vector3D.h"
#include "atlas/interpolation/element/SphericalPolygon3D.h"
#include "eckit/types/FloatCompare.h"
#include "tests/AtlasTestEnvironment.h"

using atlas::interpolation::Vector3D;
using atlas::interpolation::element::SphericalPolygon3D;
using Eigen::Vector4d;

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

const double relative_error           = 0.0001;
static const double parametricEpsilon = 1e-15;

const double recipRoot2 = 1 / sqrt(2);
const double recipRoot3 = 1 / sqrt(3);


CASE("test_spherical_polygon_triag") {
    const std::array<const Vector3D, 3> testTriangleVertices = {Vector3D{0, 1, 0}, Vector3D{recipRoot2, 0, recipRoot2},
                                                                Vector3D{1, 0, 0}};
    SphericalPolygon3D<3> testTriangle(testTriangleVertices);
    double testTriangleArea = testTriangle.area();

    SECTION("test great circle normals") {
        std::array<Vector3D, 3> expectedGreatCircleNormals = {Vector3D{recipRoot2, 0, -1 * recipRoot2},
                                                              Vector3D{0, recipRoot2, 0}, Vector3D{0, 0, 1}};
        for (size_t i = 0; i < 3; ++i) {
            EXPECT(expectedGreatCircleNormals[i].isApprox(testTriangle.getGreatCircleNormals()[i], relative_error));
        }
    }

    SECTION("test shoelace formula for area") {  // same as usual angle for triangle
        double expectedArea = 0.5210;
        Log::info() << testTriangleArea << " " << expectedArea << std::endl;
        EXPECT(eckit::types::is_approximately_equal(testTriangleArea, expectedArea, relative_error));
    }

    SECTION("test intersection and weight computation") {
        const double edgeEpsilon      = parametricEpsilon * sqrt(testTriangleArea);
        const size_t numberTestPoints = 10;

        std::array<Vector3D, numberTestPoints> candidatePoints = {
            Vector3D{0, 1, 0},                                            // Vertex 0
            Vector3D{recipRoot2, 0, recipRoot2},                          // Vertex 1
            Vector3D{1, 0, 0},                                            // Vertex 2
            Vector3D{recipRoot3, recipRoot3, recipRoot3},                 // Inside
            Vector3D{recipRoot2, recipRoot2, 0},                          // On edge
            Vector3D{recipRoot2, recipRoot2, 0.000001},                   // Just inside edge
            Vector3D{0.01, 0.9999, 0.01},                                 // Just inside vertex
            Vector3D{0, recipRoot2, recipRoot2},                          // Outside
            Vector3D{-1 * recipRoot3, -1 * recipRoot3, -1 * recipRoot3},  // Opposite Side
            Vector3D{recipRoot2, recipRoot2, -0.000001}                   // Just outside edge
        };

        std::array<bool, numberTestPoints> isPointInside = {1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
        size_t expectedInside  = std::accumulate(std::begin(isPointInside), std::end(isPointInside), 0);
        size_t expectedOutside = numberTestPoints - expectedInside;

        std::array<Vector3D, numberTestPoints> candidateWeights = {
            Vector3D{1, 0, 0},            // Vertex 0
            Vector3D{0, 1, 0},            // Vertex 1
            Vector3D{0, 0, 1},            // Vertex 2
            Vector3D{0.5774, 0.8165, 0},  // Inside
            Vector3D{0.7071, 0, 0.7071},  // On edge
            Vector3D{0.7071, 0, 0.7071},  // Just inside edge
            Vector3D{0.9999, 0.0141, 0},  // Just inside vertex
            Vector3D{0, 0, 0},            // Outside
            Vector3D{0, 0, 0},            // Opposite Side
            Vector3D{0, 0, 0}             // Just outside edge
        };

        size_t pointsInside  = 0;
        size_t pointsOutside = 0;

        for (size_t i = 0; i < numberTestPoints; ++i) {
            std::optional<std::array<double, 3>> polygonWeights =
                testTriangle.computeWeights(candidatePoints[i], edgeEpsilon);
            if (!polygonWeights) {
                EXPECT((isPointInside[i] == 0));
                pointsOutside += 1;
            }
            else {
                EXPECT((isPointInside[i] == 1));
                pointsInside += 1;
                Vector3D polygonWeightsVector((*polygonWeights)[0], (*polygonWeights)[1], (*polygonWeights)[2]);
                EXPECT(polygonWeightsVector.isApprox(candidateWeights[i], relative_error));
            }
        }
        Log::info() << "Points in/out: " << pointsInside << "/" << pointsOutside << std::endl;
        EXPECT(pointsOutside == expectedOutside);
        EXPECT(pointsInside == expectedInside);
    }
}

CASE("test_spherical_polygon_planar_quad") {
    const std::array<const Vector3D, 4> testPlanarQuadVertices = {
        Vector3D{recipRoot3, recipRoot3, recipRoot3}, Vector3D{recipRoot3, -1 * recipRoot3, recipRoot3},
        Vector3D{recipRoot3, -1 * recipRoot3, -1 * recipRoot3}, Vector3D{recipRoot3, recipRoot3, -1 * recipRoot3}};
    SphericalPolygon3D<4> testPlanarQuad(testPlanarQuadVertices);
    double testPlanarQuadArea = testPlanarQuad.area();

    SECTION("test great circle normals") {
        double twoThirds = 2. / 3;

        std::array<Vector3D, 4> expectedGreatCircleNormals = {
            Vector3D{twoThirds, 0, -1 * twoThirds}, Vector3D{twoThirds, twoThirds, 0},
            Vector3D{twoThirds, 0, twoThirds}, Vector3D{twoThirds, -1 * twoThirds, 0}};
        std::array<Vector3D, 4> greatCircleNormalsLocal = testPlanarQuad.getGreatCircleNormals();

        for (size_t i = 0; i < 4; ++i) {
            EXPECT(expectedGreatCircleNormals[i].isApprox(greatCircleNormalsLocal[i], relative_error));
        }
    }

    SECTION("test shoelace formula for area") {  // should be same as area of flat quad
        double expectedArea = 1.3333;
        Log::info() << testPlanarQuadArea << " " << expectedArea << std::endl;
        EXPECT(eckit::types::is_approximately_equal(testPlanarQuadArea, expectedArea, relative_error));
    }

    SECTION("test intersection and weight computation") {
        const double edgeEpsilon      = parametricEpsilon * sqrt(testPlanarQuadArea);
        const size_t numberTestPoints = 11;

        std::array<Vector3D, numberTestPoints> candidatePoints = {
            Vector3D{recipRoot3, recipRoot3, recipRoot3},            // Vertex 0
            Vector3D{recipRoot3, -1 * recipRoot3, recipRoot3},       // Vertex 1
            Vector3D{recipRoot3, -1 * recipRoot3, -1 * recipRoot3},  // Vertex 2
            Vector3D{recipRoot3, recipRoot3, -1 * recipRoot3},       // Vertex 3
            Vector3D{1, 0, 0},                                       // Inside
            Vector3D{recipRoot2, recipRoot2, 0},                     // On edge
            Vector3D{0.7072, 0.7070, 0},                             // Just inside edge
            Vector3D{0.5774, 0.5773, 0.5773},                        // Just inside vertex
            Vector3D{0, 1, 0},                                       // Outside
            Vector3D{-1, 0, 0},                                      // Opposite side
            Vector3D{0.7070, 0.7072, 0}                              // Just outside edge
        };

        std::array<bool, numberTestPoints> isPointInside = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
        size_t expectedInside  = std::accumulate(std::begin(isPointInside), std::end(isPointInside), 0);
        size_t expectedOutside = numberTestPoints - expectedInside;

        std::array<Vector4d, numberTestPoints> candidateWeights = {
            Vector4d{1, 0, 0, 0},
            Vector4d{0, 1, 0, 0},
            Vector4d{0, 0, 1, 0},
            Vector4d{0, 0, 0, 1},
            Vector4d{0.4330, 0.4330, 0.4330, 0.4330},
            Vector4d{0.6124, 0, 0, 0.6124},
            Vector4d{0.61237, 0.00008, 0.00008, 0.61237},
            Vector4d{1, 0, 0, 0},
            Vector4d{0, 0, 0, 0},
            Vector4d{0, 0, 0, 0},
            Vector4d{0, 0, 0, 0}  //
        };

        size_t pointsInside  = 0;
        size_t pointsOutside = 0;
        for (size_t i = 0; i < numberTestPoints; ++i) {
            std::optional<std::array<double, 4>> polygonWeights =
                testPlanarQuad.computeWeights(candidatePoints[i], edgeEpsilon);

            if (!polygonWeights) {
                EXPECT((isPointInside[i] == 0));
                pointsOutside += 1;
            }
            else {
                EXPECT((isPointInside[i] == 1));
                pointsInside += 1;
                Vector4d polygonWeightsVector((*polygonWeights)[0], (*polygonWeights)[1], (*polygonWeights)[2],
                                              (*polygonWeights)[3]);
                EXPECT(polygonWeightsVector.isApprox(candidateWeights[i], relative_error));
            }
        }
        Log::info() << "Points in/out: " << pointsInside << "/" << pointsOutside << std::endl;
        EXPECT(pointsOutside == expectedOutside);
        EXPECT(pointsInside == expectedInside);
    }
}

CASE("test_spherical_polygon_nonplanar_quad") {
    const std::array<const Vector3D, 4> testQuadVertices = {Vector3D{1, 0, 0}, Vector3D{recipRoot2, recipRoot2, 0},
                                                            Vector3D{recipRoot3, recipRoot3, recipRoot3},
                                                            Vector3D{recipRoot2, 0, recipRoot2}};
    SphericalPolygon3D<4> testQuad(testQuadVertices);
    double testQuadArea           = testQuad.area();
    const size_t numberTestPoints = 11;

    std::array<Vector3D, numberTestPoints> candidatePoints = {
        Vector3D{1, 0, 0},                             // Vertex 0
        Vector3D{recipRoot2, recipRoot2, 0},           // Vertex 1
        Vector3D{recipRoot3, recipRoot3, recipRoot3},  // Vertex2
        Vector3D{recipRoot2, 0, recipRoot2},           // Vertex 3
        Vector3D{0.8246, 0.4, 0.4},                    // Inside
        Vector3D{0.5 * sqrt(3), 0.5, 0},               // On edge
        Vector3D{0.5 * sqrt(3), 0.5, 0.000001},        // Just inside edge
        Vector3D{0.5918, 0.57, 0.57},                  // Just inside vertex
        Vector3D{recipRoot2, -1 * recipRoot2, 0},      // Outside
        Vector3D{-0.8246, -0.4, -0.4},                 // Opposite side
        Vector3D{0.5 * sqrt(3), 0.5, -0.000001},       // Just outside edge
    };

    std::array<bool, numberTestPoints> isPointInside = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};

    const size_t expectedInside  = 8;
    const size_t expectedOutside = numberTestPoints - expectedInside;

    std::array<Vector4d, expectedInside> polygonWeightsVector;

    SECTION("test great circle normals") {
        double recipRoot6 = recipRoot2 * recipRoot3;

        std::array<Vector3D, 4> expectedGreatCircleNormals = {
            Vector3D{0, 0, recipRoot2}, Vector3D{recipRoot6, -1 * recipRoot6, 0},
            Vector3D{recipRoot6, 0, -1 * recipRoot6}, Vector3D{0, recipRoot2, 0}};
        std::array<Vector3D, 4> greatCircleNormalsLocal = testQuad.getGreatCircleNormals();

        for (size_t i = 0; i < 4; ++i) {
            EXPECT(expectedGreatCircleNormals[i].isApprox(greatCircleNormalsLocal[i], relative_error));
        }
    }

    SECTION("test shoelace formula for area") {  // shadow area as quad not flat
        double expectedArea = 0.4597;
        Log::info() << testQuadArea << " " << expectedArea << std::endl;
        EXPECT(eckit::types::is_approximately_equal(testQuadArea, expectedArea, relative_error));
    }

    SECTION("test intersection and weight computation") {
        const double edgeEpsilon = parametricEpsilon * sqrt(testQuadArea);

        std::array<Vector4d, numberTestPoints> candidateWeights = {
            Vector4d{1, 0, 0, 0},
            Vector4d{0, 1, 0, 0},
            Vector4d{0, 0, 1, 0},
            Vector4d{0, 0, 0, 1},
            Vector4d{0.23847, 0.26325, 0.37043, 0.26325},
            Vector4d{0.3660, 0.7071, 0, 0},
            Vector4d{0.3660, 0.7071, 0, 0},
            Vector4d{0.0075, 0.0203, 0.9624, 0.0203},
            Vector4d{0, 0, 0, 0},
            Vector4d{0, 0, 0, 0},
            Vector4d{0, 0, 0, 0}  //
        };

        size_t pointsInside  = 0;
        size_t pointsOutside = 0;

        for (size_t i = 0; i < numberTestPoints; ++i) {
            std::optional<std::array<double, 4>> polygonWeights =
                testQuad.computeWeights(candidatePoints[i], edgeEpsilon);
            if (!polygonWeights) {
                EXPECT((isPointInside[i] == 0));
                pointsOutside += 1;
            }
            else {
                EXPECT((isPointInside[i] == 1));
                polygonWeightsVector[pointsInside] =
                    Vector4d((*polygonWeights)[0], (*polygonWeights)[1], (*polygonWeights)[2], (*polygonWeights)[3]);
                EXPECT(polygonWeightsVector[pointsInside].isApprox(candidateWeights[i], relative_error));
                pointsInside += 1;
            }
        }
        Log::info() << "Points in/out: " << pointsInside << "/" << pointsOutside << std::endl;
        EXPECT(pointsOutside == expectedOutside);
        EXPECT(pointsInside == expectedInside);
    }

    SECTION("test rotational consistency") {
        const std::array<const Vector3D, 4> testQuadRotatedVertices = {
            Vector3D{recipRoot2, 0, recipRoot2}, Vector3D{1, 0, 0}, Vector3D{recipRoot2, recipRoot2, 0},
            Vector3D{recipRoot3, recipRoot3, recipRoot3}};

        SphericalPolygon3D<4> testQuadRotated(testQuadRotatedVertices);
        double testQuadRotatedArea = testQuadRotated.area();
        const double edgeEpsilon   = parametricEpsilon * sqrt(testQuadRotatedArea);
        size_t pointsInsideRotated = 0;

        for (size_t i = 0; i < numberTestPoints; ++i) {
            std::optional<std::array<double, 4>> polygonWeightsRotated =
                testQuad.computeWeights(candidatePoints[i], edgeEpsilon);
            if (!polygonWeightsRotated) {
                EXPECT((isPointInside[i] == 0));
            }
            else {
                EXPECT((isPointInside[i] == 1));
                Vector4d polygonWeightsVectorRotated((*polygonWeightsRotated)[0], (*polygonWeightsRotated)[1],
                                                     (*polygonWeightsRotated)[2], (*polygonWeightsRotated)[3]);
                EXPECT_EQ(polygonWeightsVector[pointsInsideRotated], polygonWeightsVectorRotated);
                pointsInsideRotated += 1;
            }
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
