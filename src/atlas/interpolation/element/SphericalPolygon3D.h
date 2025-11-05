/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <iosfwd>
#include <limits>
#include <optional>

#include "atlas/interpolation/Vector3D.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace interpolation {
namespace method {
struct Ray;
}
namespace element {

//----------------------------------------------------------------------------------------------------------------------

// Spherical polygon structure
// Implements spherical mean value coordinates from
// M. Floater, “Generalized barycentric coordinates and applications” Acta Numerica, p. 001, 2016.


template <size_t numVertices>
class SphericalPolygon3D {
public:
    // Compute normals to great circles described by vertex pairs v_i x v_i+1 on construction
    SphericalPolygon3D(const std::array<const Vector3D, numVertices>& listVertices): polygonVertices(listVertices) {
        for (size_t normalsIndex = 0; normalsIndex < numVertices; ++normalsIndex) {
            greatCircleNormals[normalsIndex] =
                polygonVertices[normalsIndex].cross(polygonVertices[(normalsIndex + 1) % numVertices]);
        };
        radius = polygonVertices[0].norm();
        // TODO should probably make sure all vertices have same radius (to tolerance) and check point also
    };

    // Compute shadow area of polygon using shoelace formula
    double area() const {
        Vector3D prePolygonArea{0, 0, 0};
        for (const Vector3D& normal : greatCircleNormals) {
            prePolygonArea += normal;
        }
        return 0.5 * prePolygonArea.norm();
    }

    // Print vertices of polygon
    void print(std::ostream& stream) const {
        stream << "SphericalPolygon3D[";
        for (size_t vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) {
            stream << "vertex " << vertexIndex << ": " << polygonVertices[vertexIndex] << ", ";
        }
        stream << "]";
    }

    // Insertion operator overload
    friend std::ostream& operator<<(std::ostream& stream, const SphericalPolygon3D<numVertices>& p) {
        p.print(stream);
        return stream;
    }

    // Get vertex at index vertexIndex
    const Vector3D& p(size_t vertexIndex) {
        if (vertexIndex < numVertices) {
            return polygonVertices[vertexIndex];
        };
        throw_OutOfRange("SphericalPolygon3D::p(vertexIndex)", vertexIndex, numVertices, Here());
    }

    std::array<Vector3D, numVertices> getGreatCircleNormals() { return greatCircleNormals; }

    // Compute vertex weights in Spherical Mean Value Coordinates
    std::optional<std::array<double, numVertices>> computeWeights(
        const Vector3D& candidatePoint, double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon()) const {
        // Check if point is in polygon
        std::array<double, numVertices> greatCircleProducts = {0};

        for (size_t normalIndex = 0; normalIndex < numVertices; ++normalIndex) {
            greatCircleProducts[normalIndex] = greatCircleNormals[normalIndex].dot(candidatePoint);
            if (candidatePoint.dot(greatCircleNormals[normalIndex]) < -1 * radius * radius * radius * edgeEpsilon) {
                // multiplying by r^3 same as normalising vectors on LHS
                return {};
            }
        }

        // Compute the weights in SMV coordinates
        std::array<double, numVertices> weights = {0};

        for (size_t vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) {
            if ((candidatePoint - polygonVertices[vertexIndex]).norm() < edgeEpsilon) {
                weights[vertexIndex] = 1;
                return weights;
            }
        }

        std::array<Vector3D, numVertices> spokeNormals;
        std::array<double, numVertices> spokeNorms;
        double denominator = 0.;

        for (size_t vertexIndex = 0; vertexIndex < numVertices; ++vertexIndex) {
            spokeNormals[vertexIndex] = polygonVertices[vertexIndex].cross(candidatePoint);
            spokeNorms[vertexIndex]   = spokeNormals[vertexIndex].norm();
        }

        for (size_t weightIndex = 0; weightIndex < numVertices; ++weightIndex) {
            size_t previousWeightIndex = (numVertices + weightIndex - 1) % numVertices;
            size_t nextWeightIndex     = (weightIndex + 1) % numVertices;
            double product             = 1.;

            for (size_t index = 0; index < numVertices; ++index) {
                if ((index != weightIndex) && (index != previousWeightIndex)) {
                    product *= greatCircleProducts[index];
                }
            }

            double leftSideAngle =
                spokeNorms[nextWeightIndex] -
                spokeNormals[nextWeightIndex].dot(spokeNormals[weightIndex] / spokeNorms[weightIndex]);
            double rightSideAngle =
                spokeNorms[previousWeightIndex] -
                spokeNormals[previousWeightIndex].dot(spokeNormals[weightIndex] / spokeNorms[weightIndex]);

            weights[weightIndex] = (greatCircleProducts[previousWeightIndex] * (leftSideAngle) +
                                    greatCircleProducts[weightIndex] * (rightSideAngle)) *
                                   product;
            denominator += weights[weightIndex] * polygonVertices[weightIndex].dot(candidatePoint);
        }

        for (double& weight : weights) {
            weight *= (radius * radius);
            weight /= denominator;
        }

        return weights;
    }

private:
    // These are set on construction and do not change after
    std::array<Vector3D, numVertices> greatCircleNormals;
    const std::array<const Vector3D, numVertices>& polygonVertices;
    double radius = 1;
};


//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
