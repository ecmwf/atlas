/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <map>

#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

// #include "tests/TestMeshes.h"
#include "atlas/mpi/mpi.h"
#include "atlas/util/Debug.h"
// #include "atlas/Mesh.h"
#include "atlas/array.h"
#include "atlas/util/IndexView.h"
// #include "atlas/actions/BuildParallelFields.h"
// #include "atlas/actions/BuildPeriodicBoundaries.h"
// #include "atlas/actions/BuildHalo.h"
// #include "atlas/actions/BuildEdges.h"
// #include "atlas/actions/BuildDualMesh.h"

#include "eckit/maths/Matrix.h"
#include "tests/AtlasTestEnvironment.h"

using namespace eckit::maths;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

template <int NDIM>
class Mononomial {
public:
    Mononomial() {}

    Mononomial( double coeff_, int power_[NDIM] ) {
        coeff = coeff_;
        std::copy( power_, power_ + NDIM, power );
    }

    Mononomial( double coeff_, int p1_ ) {
        coeff    = coeff_;
        power[0] = p1_;
    }

    Mononomial( double coeff_, int p1_, int p2_ ) {
        coeff    = coeff_;
        power[0] = p1_;
        power[1] = p2_;
    }

    Mononomial( double coeff_, int p1_, int p2_, int p3_ ) {
        coeff    = coeff_;
        power[0] = p1_;
        power[1] = p2_;
        power[2] = p3_;
    }

    double eval( const double point[NDIM] ) const {
        double val = coeff;
        for ( size_t d = 0; d < NDIM; ++d ) {
            val *= std::pow( point[d], power[d] );
        }
        return val;
    }

    double eval( const double& x ) const { return coeff * std::pow( x, power[0] ); }

    double eval( const double& x, const double& y ) const {
        return coeff * std::pow( x, power[0] ) * std::pow( y, power[1] );
    }

    double eval( const double& x, const double& y, const double& z ) const {
        return coeff * std::pow( x, power[0] ) * std::pow( y, power[1] ) * std::pow( z, power[2] );
    }

    Mononomial deriv( size_t dim ) const {
        double c = coeff * power[dim];
        int p[NDIM];
        std::copy( power, power + NDIM, p );
        p[dim] -= 1;
        return Mononomial( c, p );
    }

    double coeff;
    int power[NDIM];
};

template <size_t NDIM>
class Polynomial {
public:
    typedef Mononomial<NDIM> monomial_type;

public:
    Polynomial& add( const monomial_type& mononomial ) {
        if ( mononomial.coeff != 0 )
            mononomials_.push_back( mononomial );
        return *this;
    }

    double eval( const double point[NDIM] ) const {
        double val = 0;
        for ( size_t n = 0; n < mononomials_.size(); ++n ) {
            val += mononomials_[n].eval( point );
        }
        return val;
    }

    double eval( const double& x ) const {
        double val = 0;
        for ( size_t n = 0; n < mononomials_.size(); ++n ) {
            val += mononomials_[n].eval( x );
        }
        return val;
    }

    double eval( const double& x, const double& y ) const {
        double val = 0;
        for ( size_t n = 0; n < mononomials_.size(); ++n ) {
            val += mononomials_[n].eval( x, y );
        }
        return val;
    }

    double eval( const double& x, const double& y, const double& z ) const {
        double val = 0;
        for ( size_t n = 0; n < mononomials_.size(); ++n ) {
            val += mononomials_[n].eval( x, y, z );
        }
        return val;
    }

    double operator()( const double point[NDIM] ) const { return eval( point ); }

    double operator()( const double& x ) const { return eval( x ); }

    double operator()( const double& x, const double& y ) const { return eval( x, y ); }

    double operator()( const double& x, const double& y, const double& z ) const { return eval( x, y, z ); }

    Polynomial deriv( size_t dim ) const {
        Polynomial p;
        for ( size_t i = 0; i < mononomials_.size(); ++i )
            p.add( mononomials_[i].deriv( dim ) );
        return p;
    }

    std::vector<Polynomial> grad() const {
        std::vector<Polynomial> g;
        for ( size_t d = 0; d < NDIM; ++d )
            g.push_back( deriv( d ) );
        return g;
    }

    static Polynomial div( const std::vector<Polynomial>& pvec ) {
        Polynomial p;
        for ( size_t d = 0; d < NDIM; ++d )
            p += pvec[d].deriv( d );
        return p;
    }

    static Polynomial curl_z( const std::vector<Polynomial>& pvec ) { return pvec[1].deriv( 0 ) - pvec[0].deriv( 1 ); }

    static std::vector<Polynomial> curl( const std::vector<Polynomial>& pvec ) {
        std::vector<Polynomial> p( 3 );
        if ( NDIM == 2 ) {
            p[2] = pvec[1].deriv( 0 ) - pvec[0].deriv( 1 );
        }
        if ( NDIM == 3 ) {
            p[0] = pvec[2].deriv( 1 ) - pvec[1].deriv( 2 );
            p[1] = pvec[0].deriv( 2 ) - pvec[2].deriv( 0 );
            p[2] = pvec[1].deriv( 0 ) - pvec[0].deriv( 1 );
        }
        return p;
    }

    Polynomial operator+( const Polynomial& p2 ) {
        Polynomial p = *this;
        p += p2;
        return p;
    }

    Polynomial operator-( const Polynomial& p2 ) {
        Polynomial p = *this;
        p -= p2;
        return p;
    }

    Polynomial& operator+=( const Polynomial& other ) { return addition( other, +1. ); }

    Polynomial& operator-=( const Polynomial& other ) { return addition( other, -1. ); }

private:
    Polynomial& addition( const Polynomial& other, double sign = 1. ) {
        std::vector<bool> matched( other.mononomials_.size(), false );
        for ( size_t i = 0; i < mononomials_.size(); ++i ) {
            for ( size_t j = 0; j < other.mononomials_.size(); ++j ) {
                if ( !matched[j] ) {
                    bool match = true;
                    for ( size_t d = 0; d < NDIM; ++d ) {
                        if ( mononomials_[i].power[d] != other.mononomials_[j].power[d] ) {
                            match = false;
                            break;
                        }
                    }
                    if ( match ) {
                        matched[j] = true;
                        mononomials_[i].coeff += sign * other.mononomials_[j].coeff;
                    }
                }
            }
        }
        for ( size_t j = 0; j < other.mononomials_.size(); ++j ) {
            if ( !matched[j] ) {
                add( other.mononomials_[j] );
            }
        }
        return *this;
    }

private:
    std::vector<monomial_type> mononomials_;
};

class ShapeFunction : public util::Object {
public:
    ShapeFunction() {}
    virtual ~ShapeFunction() {}
};

class ElementType : public util::Object {
public:
    typedef util::ObjectHandle<ElementType> Ptr;
    typedef std::vector<ElementType::Ptr> Vector;

public:
    size_t N() const { return N_; }

    bool parametric() const { return parametric_; }

    const double* nodes() const { return nodes_.data(); }

    const ShapeFunction& shape_function() const { return *shape_function_; }

protected:
    bool parametric_;

    size_t N_;

    std::vector<double> nodes_;

    ShapeFunction::Ptr shape_function_;
};

class Point : public ElementType {
public:
    Point() {
        N_              = 1;
        parametric_     = false;
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class Polygon : public ElementType {
public:
    Polygon( size_t max_nodes = 0 ) {
        parametric_ = false;
        N_          = max_nodes;
    }
};

class QuadP1 : public ElementType {
public:
    QuadP1() {
        parametric_      = true;
        N_               = 4;
        int nodes_init[] = {-1., -1., 1., -1., 1., 1., -1., 1.};
        nodes_.assign( nodes_init, nodes_init + N_ );
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class TriagP1 : public ElementType {
public:
    TriagP1() {
        parametric_      = true;
        N_               = 3;
        int nodes_init[] = {0., 0., 1., 0., 0., 1.};
        nodes_.assign( nodes_init, nodes_init + N_ );
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class LineP0 : public ElementType {
public:
    LineP0() {
        parametric_         = true;
        N_                  = 1;
        double nodes_init[] = {0.};
        nodes_.assign( nodes_init, nodes_init + N_ );
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class LineP1 : public ElementType {
public:
    LineP1() {
        parametric_         = true;
        N_                  = 2;
        double nodes_init[] = {-1., 1.};
        nodes_.assign( nodes_init, nodes_init + N_ );
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class LineP2 : public ElementType {
public:
    LineP2() {
        parametric_         = true;
        N_                  = 3;
        double nodes_init[] = {-1., 1., 0.};
        nodes_.assign( nodes_init, nodes_init + N_ );
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class Structured1D : public ElementType {
public:
    Structured1D( int N ) {
        parametric_     = true;
        N_              = N;
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
};

class Nodes : public util::Object {
public:
    Nodes() {
        npts_   = 0;
        nlev_   = 0;
        nlon_   = 0;
        nlat_   = 0;
        nlev_   = 0;
        nproma_ = 0;
        nblk_   = 0;
    }
    size_t npts() const { return npts_; }
    size_t nlev() const { return nlev_; }
    size_t nproma() const { return nproma_; }
    size_t nblk() const { return nblk_; }
    size_t nlon() const { return nlon_; }
    size_t nlat() const { return nlat_; }

private:
    size_t npts_;
    size_t nlev_;
    size_t nproma_;
    size_t nblk_;
    size_t nlon_;
    size_t nlat_;
};

enum IndexType
{
    IDX_NOTUSED = -100,
    IDX_NODE    = -1,
    IDX_LEVEL   = -2,
    IDX_BLK     = -3,
    IDX_NPROMA  = -4,
    IDX_VAR     = -5
};

class NewFunctionSpace {
public:
    NewFunctionSpace() : nproma_( 0 ), nb_nodes_( 0 ), nblk_( 0 ) {}

    /// @brief Number of element types
    size_t nb_element_types() const { return nelem_per_type_.size(); }

    /// @brief Element type at index
    const ElementType& element_type( size_t idx ) const {
        ATLAS_ASSERT( idx < element_types_.size() );
        return *element_types_[idx];
    }

    /// @brief Number of elements for element type at given index
    size_t nb_elements_in_element_type( size_t idx ) const {
        ATLAS_ASSERT( idx < nelem_per_type_.size() );
        return nelem_per_type_[idx];
    }

    /// @brief Number of elements for element type at given index
    size_t nb_elements() const { return nb_elements_; }

    /// @brief Number of owned elements for element type at given index
    size_t nb_owned_elements() const { return nb_owned_elements_; }

    /// @brief Number of nodes in this function space
    size_t nb_nodes() const { return nb_nodes_; }

    /// @brief Number of nodes in this function space
    size_t nb_owned_nodes() const { return nb_owned_nodes_; }

    /// @brief Number of levels in this function space
    size_t nb_levels() const { return nlev_; }

    /// @brief Number of blocks (tiles) in this function space
    size_t nblk() const { return nblk_; }

    /// @brief Number of nodes within one block in this function space
    size_t nproma() const { return nproma_; }

    /// @brief Maximum number of points per element type
    size_t N_max() const { return N_max_; }

    /// @brief Minimum number of points per element type
    size_t N_min() const { return N_min_; }

    /// @brief Set the total number of nodes in FunctionSpace
    void set_nb_nodes( int nb_nodes ) {
        nb_nodes_ = nb_nodes;
        if ( nproma_ == 0 )
            nproma_ = 1;
        ATLAS_ASSERT( nb_nodes_ % nproma_ == 0 );
        nblk_ = nb_nodes_ / nproma_;
    }

    /// @brief Set the total number of nodes in FunctionSpace
    void set_nproma( int nproma ) {
        nproma_ = nproma;
        nblk_   = 0;
        if ( nb_nodes_ != 0 ) {
            ATLAS_ASSERT( nb_nodes_ % nproma_ == 0 );
            nblk_ = nb_nodes_ / nproma_;
        }
    }

    /// @brief Set the number of levels in FunctionSpace
    void set_nb_levels( int nb_levels ) { nlev_ = nb_levels; }

    /// @brief Add elements
    const ElementType& add_elements( const std::string& element_type, int nelem ) {
        if ( element_type == "Polygon" )
            element_types_.push_back( ElementType::Ptr( new Polygon() ) );
        if ( element_type == "QuadP1" )
            element_types_.push_back( ElementType::Ptr( new QuadP1 ) );
        if ( element_type == "TriagP1" )
            element_types_.push_back( ElementType::Ptr( new TriagP1 ) );
        if ( element_type == "LineP0" )
            element_types_.push_back( ElementType::Ptr( new LineP0 ) );
        if ( element_type == "LineP1" )
            element_types_.push_back( ElementType::Ptr( new LineP1 ) );
        if ( element_type == "LineP2" )
            element_types_.push_back( ElementType::Ptr( new LineP2 ) );
        index_[element_type] = nelem_per_type_.size();
        nelem_per_type_.push_back( nelem );

        const ElementType& elem = *element_types_.back();

        if ( nelem_per_type_.size() == 1 ) {
            N_max_       = elem.N();
            N_min_       = elem.N();
            nb_elements_ = nelem;
        }
        else {
            N_max_ = std::max( N_max_, elem.N() );
            N_min_ = std::min( N_min_, elem.N() );
            nb_elements_ += nelem;
        }
        return elem;
    }

    /// @brief Set the nodes that define the element
    void set_nodes( const Nodes& nodes ) { nodes_ = &nodes; }

    template <typename DATA_TYPE>
    atlas::array::ArrayT<DATA_TYPE> create_field( int idx1 = IDX_NOTUSED, int idx2 = IDX_NOTUSED,
                                                  int idx3 = IDX_NOTUSED, int idx4 = IDX_NOTUSED,
                                                  int idx5 = IDX_NOTUSED, int idx6 = IDX_NOTUSED,
                                                  int idx7 = IDX_NOTUSED ) {
        atlas::array::ArrayShape shape;
        shape.reserve( 7 );
        if ( idx1 != IDX_NOTUSED )
            shape.push_back( range( idx1 ) );
        if ( idx2 != IDX_NOTUSED )
            shape.push_back( range( idx2 ) );
        if ( idx3 != IDX_NOTUSED )
            shape.push_back( range( idx3 ) );
        if ( idx4 != IDX_NOTUSED )
            shape.push_back( range( idx4 ) );
        if ( idx5 != IDX_NOTUSED )
            shape.push_back( range( idx5 ) );
        if ( idx6 != IDX_NOTUSED )
            shape.push_back( range( idx6 ) );
        if ( idx7 != IDX_NOTUSED )
            shape.push_back( range( idx7 ) );
        atlas::array::ArrayT<DATA_TYPE> field( shape );
        return field;
    }

private:
    int range( int idx_type ) const {
        switch ( idx_type ) {
            case IDX_NODE:
                return nb_nodes();
            case IDX_LEVEL:
                return nb_levels();
            case IDX_BLK:
                return nblk();
            case IDX_NPROMA:
                return nproma();
            default:
                if ( idx_type >= 0 )
                    return idx_type;
        }
        throw_Exception( "idx_type not recognized" );
        return 0;
    }

private:
    /// @brief Lookup of element_type index by name
    std::map<std::string, int> index_;

    /// @brief Reference to nodes
    Nodes const* nodes_;

    /// @brief Element types
    ElementType::Vector element_types_;

    /// @brief Number of elements per element type
    std::vector<int> nelem_per_type_;

    /// @brief Total number of elements
    size_t nb_elements_;

    /// @brief Maximum number of points per element type
    size_t N_max_;

    /// @brief Minimum number of points per element type
    size_t N_min_;

    /// @brief Total number of nodes
    size_t nb_nodes_;

    /// @brief Number of vertical levels
    size_t nlev_;

    /// @brief Number of nodes in one blk
    size_t nproma_;

    /// @brief Number of blocks in domain
    size_t nblk_;

    size_t nb_owned_nodes_;
    size_t nb_owned_elements_;
};

class Column : public ElementType  // Really is a Elementtype in 3D
{
public:
    Column( const ElementType::Ptr& horizontal, const ElementType::Ptr& vertical ) :
        horizontal_( horizontal ), vertical_( vertical ) {
        shape_function_ = ShapeFunction::Ptr( new ShapeFunction );
    }
    size_t N() const { return npts() * nlev(); }
    size_t npts() const { return horizontal_->N(); }
    size_t nlev() const { return vertical_->N(); }

protected:
    ElementType::Ptr horizontal_;
    ElementType::Ptr vertical_;
};

template <int NDIM>
std::vector<Polynomial<NDIM>> polynomial_basis( int order, double* points, int npts ) {
    std::vector<Polynomial<NDIM>> basis( npts );

    Matrix<double> vandermonde( npts, npts );
    Matrix<double> coefficients( npts, npts );
    Matrix<int> powers( npts, NDIM );
    Matrix<double>::Proxy pts( points, npts, NDIM );
    if ( NDIM == 1 ) {
        size_t n( 0 );
        for ( size_t k = 0; k <= order; ++k ) {
            powers( n, 0 ) = k;
            ++n;
        }
        for ( n = 0; n < npts; ++n ) {
            double x = pts( n, 0 );
            for ( size_t p = 0; p < npts; ++p ) {
                double x_p          = ( powers( p, 0 ) == 0. ? 1. : std::pow( x, powers( p, 0 ) ) );
                vandermonde( n, p ) = x_p;
            }
        }
    }
    else if ( NDIM == 2 ) {
        size_t n( 0 );
        for ( size_t l = 0; l <= order; ++l ) {
            for ( size_t k = 0; k <= order; ++k ) {
                powers( n, 0 ) = k;
                powers( n, 1 ) = l;
                ++n;
            }
        }
        for ( n = 0; n < npts; ++n ) {
            double x = pts( n, 0 );
            double y = pts( n, 1 );
            for ( size_t p = 0; p < npts; ++p ) {
                double x_p          = ( powers( p, 0 ) == 0. ? 1. : std::pow( x, powers( p, 0 ) ) );
                double y_p          = ( powers( p, 1 ) == 0. ? 1. : std::pow( y, powers( p, 1 ) ) );
                vandermonde( n, p ) = x_p * y_p;
            }
        }
    }
    std::cout << "pts = \n" << pts << std::endl;
    std::cout << "powers = \n" << powers << std::endl;
    std::cout << "vandermonde = \n" << vandermonde << std::endl;
    coefficients = vandermonde.inverse().transpose();
    std::cout << "coefficients = \n" << coefficients << std::endl;
    int pwr[3];
    for ( size_t n = 0; n < npts; ++n ) {
        for ( size_t i = 0; i < npts; ++i ) {
            for ( size_t d = 0; d < NDIM; ++d )
                pwr[d] = powers( i, d );
            basis[n].add( Mononomial<NDIM>( coefficients( n, i ), pwr ) );
        }
    }
    return basis;
}

template <typename DATA_TYPE>
IndexView<DATA_TYPE, 2> make_IndexView( atlas::array::ArrayT<DATA_TYPE>& array, NewFunctionSpace& elements,
                                        int element_type_index ) {
    IndexView<DATA_TYPE, 2> view;
    size_t offset = 0;
    size_t strides[2];
    size_t shape[2];
    ATLAS_ASSERT( element_type_index < elements.nb_element_types() );
    ATLAS_ASSERT( array.shape( 1 ) == elements.N_max() );

    for ( int i = 0; i < element_type_index; ++i ) {
        offset += elements.nb_elements_in_element_type( i ) * elements.N_max();
    }

    const ElementType& elem_type = elements.element_type( element_type_index );
    strides[0]                   = array.shape( 1 );
    strides[1]                   = 1;
    shape[0]                     = elements.nb_elements_in_element_type( element_type_index );
    shape[1]                     = elem_type.N();
    return IndexView<DATA_TYPE, 2>( array.data() + offset, strides, shape );
}

class Connectivity : public IndexView<int, 2> {
    Connectivity( NewFunctionSpace& elements, int element_type_index ) {}
};

//-----------------------------------------------------------------------------

CASE( "test_functionspace" ) {
    ElementType::Ptr point( new Point );
    ElementType::Ptr quad( new QuadP1 );
    ElementType::Ptr levels( new Structured1D( 100 ) );

    Column quad_column( quad, levels );
    Column point_column( point, levels );

    // std::cout << quad->N() << std::endl;
    // std::cout << levels->N() << std::endl;
    DEBUG_VAR( quad_column.N() );
    DEBUG_VAR( quad_column.nlev() );
    DEBUG_VAR( quad_column.npts() );

    DEBUG_VAR( point_column.N() );
    DEBUG_VAR( point_column.nlev() );
    DEBUG_VAR( point_column.npts() );

    Nodes nodes;

    NewFunctionSpace fs;
    fs.set_nb_nodes( 8 );
    fs.set_nproma( 4 );
    fs.set_nb_levels( 100 );
    const ElementType& triags = fs.add_elements( "TriagP1", 2 );
    EXPECT( fs.nb_element_types() == 1 );
    EXPECT( fs.N_max() == 3 );
    EXPECT( fs.N_min() == 3 );
    EXPECT( fs.nb_elements() == 2 );
    const ElementType& quads = fs.add_elements( "QuadP1", 2 );
    EXPECT( fs.nb_element_types() == 2 );
    EXPECT( fs.N_max() == 4 );
    EXPECT( fs.N_min() == 3 );
    EXPECT( fs.nb_elements() == 4 );

    /// Allocate array for all connectivity across all elements
    atlas::array::ArrayT<int> element_node_connectivity( fs.nb_elements(), fs.N_max() );

    /// Access the data across all elements
    atlas::IndexView<int, 2> elem_connectivity( element_node_connectivity );
    // --- Triangle 1 ---
    elem_connectivity( 0, 0 ) = 1;
    elem_connectivity( 0, 1 ) = 2;
    elem_connectivity( 0, 2 ) = 3;
    elem_connectivity( 0, 3 ) = -1;
    // --- Triangle 2---
    elem_connectivity( 1, 0 ) = 11;
    elem_connectivity( 1, 1 ) = 12;
    elem_connectivity( 1, 2 ) = 13;
    elem_connectivity( 1, 3 ) = -1;
    // --- Quad 1 ---
    elem_connectivity( 2, 0 ) = 21;
    elem_connectivity( 2, 1 ) = 22;
    elem_connectivity( 2, 2 ) = 23;
    elem_connectivity( 2, 3 ) = 24;
    // --- Quad 2 ---
    elem_connectivity( 3, 0 ) = 31;
    elem_connectivity( 3, 1 ) = 32;
    elem_connectivity( 3, 2 ) = 33;
    elem_connectivity( 3, 3 ) = 34;

    /// Access the data
    atlas::IndexView<idx_t, 2> triag_node_connectivity = make_IndexView( element_node_connectivity, fs, 0 );
    EXPECT( triag_node_connectivity.shape( 0 ) == 2 );
    EXPECT( triag_node_connectivity.shape( 1 ) == 3 );
    EXPECT( triag_node_connectivity( 0, 0 ) == 1 );
    EXPECT( triag_node_connectivity( 0, 1 ) == 2 );
    EXPECT( triag_node_connectivity( 0, 2 ) == 3 );
    EXPECT( triag_node_connectivity( 1, 0 ) == 11 );
    EXPECT( triag_node_connectivity( 1, 1 ) == 12 );
    EXPECT( triag_node_connectivity( 1, 2 ) == 13 );

    BOOST_CHECK_THROW( triag_node_connectivity( 0, 3 ),
                       eckit::OutOfRange );  // should fail (OUT OF RANGE)

    atlas::IndexView<idx_t, 2> quad_node_connectivity = make_IndexView( element_node_connectivity, fs, 1 );
  EXPECT( quad_node_connectivity.shape(0 == 2 );
  EXPECT( quad_node_connectivity.shape(1 == 4 );
  EXPECT( quad_node_connectivity(0,0 == 21 );
  EXPECT( quad_node_connectivity(0,1 == 22 );
  EXPECT( quad_node_connectivity(0,2 == 23 );
  EXPECT( quad_node_connectivity(0,3 == 24 );
  EXPECT( quad_node_connectivity(1,0 == 31 );
  EXPECT( quad_node_connectivity(1,1 == 32 );
  EXPECT( quad_node_connectivity(1,2 == 33 );
  EXPECT( quad_node_connectivity(1,3 == 34 );

  NewFunctionSpace edges;
  // edges.set_nodes( nodes );
  edges.add_elements("LineP0",11);

  atlas::array::ArrayT<double> press_ifs(fs.nproma(),fs.nb_levels(),nodes.nblk());
  atlas::array::ArrayT<double> press = fs.create_field<double>(IDX_LEVEL,IDX_NODE);
  EXPECT( press.size() == fs.nb_levels()*fs.nb_nodes() );

  atlas::array::ArrayT<double> wind_uv = fs.create_field<double>(2,IDX_LEVEL,IDX_NODE);
  EXPECT( wind_uv.size() == 2*fs.nb_levels()*fs.nb_nodes() );

  EXPECT( fs.nproma() == 4);
  EXPECT( fs.nblk() == 2);
  EXPECT( fs.nb_nodes() == 8);
  atlas::array::ArrayT<double> wind_uv_ifs = fs.create_field<double>(IDX_NPROMA,IDX_LEVEL,2,IDX_BLK);
  EXPECT( wind_uv_ifs.size() == 2*fs.nb_levels()*fs.nb_nodes() );

  EXPECT( wind_uv_ifs.rank() == 4);
  EXPECT( wind_uv_ifs.shape(0) == fs.nproma() );
  EXPECT( wind_uv_ifs.shape(1) == fs.nb_levels() );
  EXPECT( wind_uv_ifs.shape(2) == 2 );
  EXPECT( wind_uv_ifs.shape(3) == fs.nblk() );
}

CASE( "test_polynomial" ) {
    typedef Polynomial<2> polynomial_type;
    typedef polynomial_type::monomial_type monomial_type;

    polynomial_type p;

    p.add( monomial_type( 3., 0, 0 ) );
    p.add( monomial_type( 1., 1, 0 ) );
    p.add( monomial_type( 1., 2, 0 ) );
    DEBUG_VAR( p( 2, 2 ) );

    polynomial_type dpdx = p.deriv( 0 );
    DEBUG_VAR( dpdx( 2, 2 ) );

    polynomial_type dpdy = p.deriv( 1 );
    DEBUG_VAR( dpdy( 2, 2 ) );

    std::vector<polynomial_type> grad = p.grad();

    std::vector<polynomial_type> pvec( 2, p );
    polynomial_type div = polynomial_type::div( pvec );
    DEBUG_VAR( div( 2, 2 ) );

    DEBUG_VAR( polynomial_type::curl_z( grad )( 2, 2 ) );

    typedef std::vector<Polynomial<1>> PolynomialBasis1D;
    typedef std::vector<Polynomial<2>> PolynomialBasis2D;

    Matrix<int> m;
    //  EXPECT( m.data_ == NULL );
    //  EXPECT( m.nr_ == 0 );
    //  EXPECT( m.nc_ == 0 );
    m.resize( 2, 3 );

    EXPECT( m.size() == 6 );
    EXPECT( m.rows() == 2 );
    EXPECT( m.cols() == 3 );

    m( 0, 0 ) = 0;
    m( 0, 1 ) = 2;
    m( 0, 2 ) = 4;
    m( 1, 0 ) = 1;
    m( 1, 1 ) = 3;
    m( 1, 2 ) = 5;

    m *= 10;

    Matrix<int> b = m + m;

    std::cout << "m = \n" << b << std::endl;

    double line_pts[]            = {-1., 1.};
    int line_npts                = 2;
    PolynomialBasis1D line_basis = polynomial_basis<1>( 1, line_pts, line_npts );
    for ( size_t n = 0; n < line_npts; ++n )
        DEBUG( n << ":  " << line_basis[n]( -1. ) << "  " << line_basis[n]( +1. ) );

    double quad_pts[]            = {-1, 1., 1., -1, -1., -1., 1., 1.};
    int quad_npts                = 4;
    PolynomialBasis2D quad_basis = polynomial_basis<2>( 1, quad_pts, quad_npts );

    for ( size_t n = 0; n < quad_npts; ++n )
        DEBUG( n << ":  " << quad_basis[n]( -1., -1. ) << "  " << quad_basis[n]( 1., -1. ) << " "
                 << quad_basis[n]( 1., 1. ) << " " << quad_basis[n]( -1., 1. ) );

    RowVector<int> vr( 3 );
    vr( 0 ) = 0;
    vr( 1 ) = 1;
    vr( 2 ) = 2;

    ColVector<int> vc( 3 );
    vc( 0 ) = 0;
    vc( 1 ) = 1;
    vc( 2 ) = 2;

    std::cout << vr << "\n" << vc << "\n\n" << vr * vc << std::endl;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
