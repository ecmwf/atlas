#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <fstream>

#include "stripack/stripack.h"

template <typename T>
inline std::ostream& operator <<(std::ostream& s, const std::vector<T>& t) {

    s << '{';
    for (size_t i = 0; i < t.size(); i++) {
        if (i != 0)
            s << ',';
        s << t[i];
    }
    s << '}';
    return s;
}
template <typename T, size_t N>
inline std::ostream& operator <<(std::ostream& s, const std::array<T,N>& t) {

    s << '{';
    for (size_t i = 0; i < t.size(); i++) {
        if (i != 0)
            s << ',';
        s << t[i];
    }
    s << '}';
    return s;
}

std::vector<std::array<double,2>> points_1();

void write_gmsh( const std::string& file, const std::vector<std::array<double,2>>& points, const std::vector<std::array<int,3>>& connectivity );

int main() {
    auto points = points_1();

    stripack::Triangulation t(points.size(),points,true);

    write_gmsh("stripack.msh",points,t.triangles());
}


class NormaliseLongitude {
public:
    // Normalise longitude between (west - eps, east - eps ) with west = 0., east = 360.
    constexpr NormaliseLongitude(): west_(-eps_), east_(360. - eps_) {}

    // Normalise longitude between ( west-eps, east-eps )  with east = west + 360
    constexpr NormaliseLongitude(double west): west_(west - eps_), east_(west + 360. - eps_) {}

    // Normalise longitude between ( west-eps, east+eps )
    constexpr NormaliseLongitude(double west, double east): west_(west - eps_), east_(east + eps_) {}

    constexpr NormaliseLongitude(const NormaliseLongitude& other): west_(other.west_), east_(other.east_) {}

    double operator()(double lon) const {
        while (lon < west_) {
            lon += 360.;
        }
        while (lon > east_) {
            lon -= 360.;
        }
        return lon;
    }

private:
    double west_;
    double east_;

public:
    static constexpr double eps_ = 1.e-11;
};

void lonlat2xyz(double lon, double lat, double& x, double& y, double& z) {
    /*
     * See https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
     * numerical conditioning for both ϕ (poles) and λ (Greenwich/Date Line).
     *
     * cos α = sqrt( 1 - sin^2 α) is better conditioned than explicit cos α, and
     * coupled with λ in [-180°, 180°[ the accuracy of the trigonometric
     * functions is the same (before converting/multiplying its angle argument
     * to radian) and explicitly chosing -180° over 180° for longitude.
     *
     * These three conditionings combined project very accurately to the sphere
     * poles and quadrants.
     */
    constexpr double degrees_to_radians = M_PI / 180.;
    constexpr NormaliseLongitude normalise_longitude(-180.);
    const double lambda_deg = normalise_longitude(lon);
    const double lambda     = degrees_to_radians * lambda_deg;
    const double phi        = degrees_to_radians * lat;

    const double sin_phi    = std::sin(phi);
    const double cos_phi    = std::sqrt(1. - sin_phi * sin_phi);
    const double sin_lambda = std::abs(lambda_deg) < 180. ? std::sin(lambda) : 0.;
    const double cos_lambda = std::abs(lambda_deg) > 90. ? std::cos(lambda) : std::sqrt(1. - sin_lambda * sin_lambda);

    x = cos_phi * cos_lambda;
    y = cos_phi * sin_lambda;
    z = sin_phi;
}

void write_gmsh(const std::string& file, const std::vector<std::array<double,2>>& nodes, const std::vector<std::array<int,3>>& elements ) {
    std::ofstream out(file.c_str());
    int gmsh_triangle_type = 2;
    out << "$MeshFormat\n";
    out << "2.2 0 " << sizeof(double) << "\n";
    out << "$EndMeshFormat\n";
    out << "$Nodes\n";
    out << nodes.size() << "\n";
    double lon, lat, x, y, z;
    for (int n = 0; n < nodes.size(); ++n) {
        int g = n+1;
        lonlat2xyz(nodes[n][0], nodes[n][1], x, y, z);
        out << g << " " << x << " " << y << " " << z << "\n";
    }
    out << "$EndNodes\n";
    out << "$Elements\n";
    out << elements.size() << "\n";
    std::string elem_info = " 2 4 1 1 1 1";
    for (int e = 0; e < elements.size(); ++e) {
        out << e+1 << elem_info;
        for (int n = 0; n < 3; ++n) {
            out << " " << elements[e][n]+1;
        }
        out << '\n';
    }
    out << "$EndElements\n";
}


std::vector<std::array<double,2>> points_1() {
    return std::vector<std::array<double,2>> {
    {180,0},
    {90,0},
    {-90,0},
    {0,90},
    {0,-90},
    {0,0},
    {18,0},
    {36,0},
    {54,0},
    {72,0},
    {108,0},
    {126,0},
    {144,0},
    {162,0},
    {-162,0},
    {-144,0},
    {-126,0},
    {-108,0},
    {-72,0},
    {-54,0},
    {-36,0},
    {-18,0},
    {0,18},
    {0,36},
    {0,54},
    {0,72},
    {180,72},
    {180,54},
    {180,36},
    {180,18},
    {180,-18},
    {180,-36},
    {180,-54},
    {180,-72},
    {0,-72},
    {0,-54},
    {0,-36},
    {0,-18},
    {90,18},
    {90,36},
    {90,54},
    {90,72},
    {-90,72},
    {-90,54},
    {-90,36},
    {-90,18},
    {-90,-18},
    {-90,-36},
    {-90,-54},
    {-90,-72},
    {90,-72},
    {90,-54},
    {90,-36},
    {90,-18},
    {123.974,-58.6741},
    {154.087,-16.9547},
    {154.212,-58.8675},
    {114.377,-41.9617},
    {125.567,-23.5133},
    {137.627,-40.8524},
    {106.162,-24.5874},
    {158.508,-38.55},
    {137.826,-72.8109},
    {142.103,-26.799},
    {138.256,-13.8871},
    {168.39,-24.3266},
    {168.954,-12.0094},
    {117.333,-12.35},
    {102.254,-11.1537},
    {120.307,59.7167},
    {107.196,26.0167},
    {144.768,28.3721},
    {150.891,60.0343},
    {164.566,25.5053},
    {116.851,14.0295},
    {124.84,28.3978},
    {157.901,42.042},
    {111.41,43.1056},
    {134.333,44.6677},
    {103.277,11.707},
    {135.358,73.2119},
    {135.349,14.2311},
    {153.48,13.386},
    {168.071,11.5344},
    {-162.99,26.3775},
    {-147.519,56.1313},
    {-122.579,27.4824},
    {-117.909,59.2376},
    {-104.052,27.3616},
    {-153.107,14.9717},
    {-110.833,41.7436},
    {-144.847,32.8534},
    {-161.546,42.1031},
    {-129.866,44.5201},
    {-133.883,72.4163},
    {-166.729,11.8907},
    {-135.755,15.2529},
    {-106.063,14.4869},
    {-119.452,11.7037},
    {-146.026,-58.6741},
    {-115.913,-16.9547},
    {-115.788,-58.8675},
    {-155.623,-41.9617},
    {-144.433,-23.5133},
    {-132.373,-40.8524},
    {-163.838,-24.5874},
    {-111.492,-38.55},
    {-132.174,-72.8109},
    {-127.897,-26.799},
    {-131.744,-13.8871},
    {-101.61,-24.3266},
    {-101.046,-12.0094},
    {-152.667,-12.35},
    {-167.746,-11.1537},
    {-14.0127,-27.2963},
    {-59.193,-57.0815},
    {-56.465,-19.5751},
    {-27.056,-59.3077},
    {-57.124,-35.9752},
    {-33.4636,-28.3914},
    {-74.8037,-46.8602},
    {-40.089,-45.1376},
    {-74.8149,-28.3136},
    {-21.3072,-42.2177},
    {-44.0778,-72.6353},
    {-19.6969,-12.8527},
    {-40.1318,-12.1601},
    {-72.691,-11.4129},
    {-56.0261,58.6741},
    {-25.9127,16.9547},
    {-25.7876,58.8675},
    {-65.6229,41.9617},
    {-54.4335,23.5133},
    {-42.373,40.8524},
    {-73.838,24.5874},
    {-21.4917,38.55},
    {-42.1744,72.8109},
    {-37.8974,26.799},
    {-41.7437,13.8871},
    {-11.6095,24.3266},
    {-11.0459,12.0094},
    {-62.667,12.35},
    {-77.7456,11.1537},
    {30.3071,59.7167},
    {17.1956,26.0167},
    {54.7676,28.3721},
    {60.8915,60.0343},
    {74.5657,25.5053},
    {26.8506,14.0295},
    {34.8398,28.3978},
    {67.9014,42.042},
    {21.41,43.1056},
    {44.3335,44.6677},
    {13.2772,11.707},
    {45.3579,73.2119},
    {45.3492,14.2311},
    {63.4799,13.386},
    {78.0714,11.5344},
    {17.01,-26.3775},
    {32.4806,-56.1313},
    {57.4213,-27.4824},
    {62.0912,-59.2376},
    {75.9483,-27.3616},
    {26.893,-14.9717},
    {69.1672,-41.7436},
    {35.1527,-32.8534},
    {18.4543,-42.1031},
    {50.1344,-44.5201},
    {46.1172,-72.4163},
    {13.2711,-11.8907},
    {44.2448,-15.2529},
    {73.9368,-14.4869},
    {60.5478,-11.7037}};
}