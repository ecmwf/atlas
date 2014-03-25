/*
 * (C) Copyright 1996-2012 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/runtime/Tool.h"

#include "eckit/grid/Grid.h"
#include "eckit/grid/Gaussian.h"

#include <math.h>
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))

using namespace std;
using namespace eckit;

//-----------------------------------------------------------------------------

namespace eckit_test {

//-----------------------------------------------------------------------------

class TestGaussian : public Tool {
public:

    TestGaussian(int argc,char **argv): Tool(argc,argv) {}

    ~TestGaussian() {}

    virtual void run();

    void test_constructor();
    void test_latitudes(const std::vector<double>& ref_data, int gaussian_number);
};

//-----------------------------------------------------------------------------

void TestGaussian::test_constructor()
{
    using namespace eckit::grid;

    BoundBox2D earth ( Point2D(-90.,0.), Point2D(90.,360.) );
    Grid* g = NULL;

    // standard case

    g = new Gaussian( 48, earth );

    ASSERT( g );
    /// @todo review this: we wrap from 0 to 360 inclusive as we do for latlon
    ASSERT( g->coordinates().size() == (48 * 2 * ( 48 * 4 + 1) ) );

    /// @todo substitute these comparisons with proper floating point comparisons
    ASSERT( g->boundingBox().bottom_left_.lat_ == -90. );
    ASSERT( g->boundingBox().bottom_left_.lon_ ==   0. );
    ASSERT( g->boundingBox().top_right_.lat_ ==  90. );
    ASSERT( g->boundingBox().top_right_.lon_ == 360. );

}


void TestGaussian::test_latitudes(const std::vector<double>& ref_data, int gaussian_number)
{
    using namespace eckit::grid;

    BoundBox2D earth ( Point2D(-90.,0.), Point2D(90.,360.) );
    Grid* g = NULL;

    // standard case

    //const int gaussian_number = 48;
    g = new Gaussian( gaussian_number, earth );
    ASSERT( g );


    // number of decimal places to compare numbers (see reference data above)
    int NDP = 6;

    ASSERT( g->coordinates().size() == (gaussian_number * 2 * ( gaussian_number * 4 + 1) ) );

    for (unsigned int i = 0; i < g->coordinates().size(); i++)
    {
        double generated_latitude = g->coordinates()[i].lat_;
        int gen = NINT(pow(10.0, NDP) * generated_latitude) ;

        // we can't be sure of the data order and don't want to enforce
        // one, so despite being slow we check whether this is a valid
        // latitude by iterating through the reference list till we find a
        // match
        bool found = false;
        for (unsigned int k = 0; k < ref_data.size(); k++)
        {
            double reference_latitude = ref_data[k];
            
            int ref = NINT(pow(10.0, NDP) * reference_latitude);
            // compare to the number of DP specified
            if (gen == ref)
            {
                found = true;
                break;
            }

        }   
        
        // if not found we have a problem
        ASSERT(found);
    }
}

//-----------------------------------------------------------------------------

void TestGaussian::run()
{
    test_constructor();



    // data from
    // http://www.ecmwf.int/publications/manuals/libraries/interpolation/n48FIS.html
    // quoted to 6 decimal places
    double g48[] = { 88.572169, 86.722531, 84.861970, 82.998942, 81.134977, 79.270559, 77.405888, 75.541061, 
                     73.676132, 71.811132, 69.946081, 68.080991, 66.215872, 64.350730, 62.485571, 60.620396, 
                     58.755209, 56.890013, 55.024808, 53.159595, 51.294377, 49.429154, 47.563926, 45.698694, 
                     43.833459, 41.968220, 40.102979, 38.237736, 36.372491, 34.507243, 32.641994, 30.776744, 
                     28.911492, 27.046239, 25.180986, 23.315731, 21.450475, 19.585219, 17.719962, 15.854704,
                     13.989446, 12.124187, 10.258928,  8.393669,  6.528409,  4.663150,  2.797890,  0.932630, 
                     -0.932630, -2.797890, -4.663150, -6.528409, -8.393669,-10.258928,-12.124187,-13.989446,
                    -15.854704,-17.719962,-19.585219,-21.450475,-23.315731,-25.180986,-27.046239,-28.911492,
                    -30.776744,-32.641994,-34.507243,-36.372491,-38.237736,-40.102979,-41.968220,-43.833459,
                    -45.698694,-47.563926,-49.429154,-51.294377,-53.159595,-55.024808,-56.890013,-58.755209,
                    -60.620396,-62.485571,-64.350730,-66.215872,-68.080991,-69.946081,-71.811132,-73.676132,
                    -75.541061,-77.405888,-79.270559,-81.134977,-82.998942,-84.861970,-86.722531,-88.572169 };


    double g80[] = {89.141519, 88.029429, 86.910771, 85.790629, 84.669924, 83.548947, 82.427818, 81.306595, 80.185310, 79.063982, 
                    77.942624, 76.821243, 75.699844, 74.578432, 73.457008, 72.335576, 71.214136, 70.092690, 68.971240, 67.849784, 
                    66.728326, 65.606864, 64.485399, 63.363932, 62.242462, 61.120991, 59.999518, 58.878044, 57.756569, 56.635092, 
                    55.513614, 54.392135, 53.270655, 52.149175, 51.027694, 49.906212, 48.784729, 47.663246, 46.541763, 45.420279, 
                    44.298794, 43.177309, 42.055824, 40.934338, 39.812852, 38.691366, 37.569880, 36.448393, 35.326906, 34.205418, 
                    33.083931, 31.962443, 30.840955, 29.719467, 28.597979, 27.476491, 26.355002, 25.233514, 24.112025, 22.990536, 
                    21.869047, 20.747558, 19.626069, 18.504580, 17.383091, 16.261601, 15.140112, 14.018622, 12.897133, 11.775643,   
                    10.654153,  9.532664,  8.411174,  7.289684,  6.168194,  5.046704,  3.925215,  2.803725,  1.682235,  0.560745,  
                    -0.560745, -1.682235, -2.803725, -3.925215, -5.046704, -6.168194, -7.289684, -8.411174, -9.532664,-10.654153,
                   -11.775643,-12.897133,-14.018622,-15.140112,-16.261601,-17.383091,-18.504580,-19.626069,-20.747558,-21.869047,
                   -22.990536,-24.112025,-25.233514,-26.355002,-27.476491,-28.597979,-29.719467,-30.840955,-31.962443,-33.083931,
                   -34.205418,-35.326906,-36.448393,-37.569880,-38.691366,-39.812852,-40.934338,-42.055824,-43.177309,-44.298794,
                   -45.420279,-46.541763,-47.663246,-48.784729,-49.906212,-51.027694,-52.149175,-53.270655,-54.392135,-55.513614,
                   -56.635092,-57.756569,-58.878044,-59.999518,-61.120991,-62.242462,-63.363932,-64.485399,-65.606864,-66.728326,
                   -67.849784,-68.971240,-70.092690,-71.214136,-72.335576,-73.457008,-74.578432,-75.699844,-76.821243,-77.942624,
                   -79.063982,-80.185310,-81.306595,-82.427818,-83.548947,-84.669924,-85.790629,-86.910771,-88.029429,-89.141519 };

    test_latitudes(std::vector<double>(g48, g48 +  96), 48);
    test_latitudes(std::vector<double>(g80, g80 + 160), 80);

}

//-----------------------------------------------------------------------------

} // namespace eckit_test

//-----------------------------------------------------------------------------

int main(int argc,char **argv)
{
    eckit_test::TestGaussian mytest(argc,argv);
    mytest.start();
    return 0;
}

