/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestFunctionSpace
#include "ecbuild/boost_test_framework.h"

#include <map>

#include "tests/TestMeshes.h"
#include "atlas/mpi/mpi.h"
#include "atlas/io/Gmsh.h"
#include "atlas/util/Debug.h"
#include "atlas/Mesh.h"
#include "atlas/util/Array.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildDualMesh.h"

#include "eckit/maths/Matrix.h"

using namespace eckit::maths;

namespace atlas {
namespace test {

  template< int NDIM >
  class Mononomial
	{
	public:
		Mononomial() {}

		Mononomial(double coeff_, int power_[NDIM] )
		{
			coeff=coeff_;
			std::copy(power_,power_+NDIM,power);
		}

		Mononomial(double coeff_, int p1_)
		{
			coeff=coeff_;
			power[0]=p1_;
		}

		Mononomial(double coeff_, int p1_, int p2_)
		{
			coeff=coeff_;
			power[0]=p1_;
			power[1]=p2_;
		}

		Mononomial(double coeff_, int p1_, int p2_, int p3_)
		{
			coeff=coeff_;
			power[0]=p1_;
			power[1]=p2_;
			power[2]=p3_;
		}

		double eval( const double point[NDIM] ) const
		{
			double val = coeff;
			for( int d=0; d<NDIM; ++d )
			{
				val *= std::pow(point[d], power[d]);
			}
			return val;
		}

		double eval( const double& x ) const
		{
			return coeff * std::pow(x,power[0]);
		}

		double eval( const double& x, const double& y ) const
		{
			return coeff * std::pow(x,power[0]) * std::pow(y,power[1]);
		}

		double eval( const double& x, const double& y, const double& z ) const
		{
			return coeff * std::pow(x,power[0]) * std::pow(y,power[1]) * std::pow(z,power[2]);
		}

		Mononomial deriv(int dim) const
		{
			double c = coeff*power[dim];
			int p[NDIM];
			std::copy(power,power+NDIM,p);
			p[dim] -= 1;
			return Mononomial( c, p );
		}

		double coeff;
		int power[NDIM];
	};

	template< int NDIM >
	class Polynomial
	{
	public:
		typedef Mononomial<NDIM> monomial_type;
	public:

		Polynomial& add( const monomial_type& mononomial )
		{
			if( mononomial.coeff != 0 )
				mononomials_.push_back(mononomial);
			return *this;
		}

		double eval( const double point[NDIM] ) const
		{
			double val = 0;
			for( int n=0; n<mononomials_.size(); ++n )
			{
				val += mononomials_[n].eval(point);
			}
			return val;
		}

		double eval( const double& x ) const
		{
			double val = 0;
			for( int n=0; n<mononomials_.size(); ++n )
			{
				val += mononomials_[n].eval(x);
			}
			return val;
		}

		double eval( const double& x, const double& y ) const
		{
			double val = 0;
			for( int n=0; n<mononomials_.size(); ++n )
			{
				val += mononomials_[n].eval(x,y);
			}
			return val;
		}

		double eval( const double& x, const double& y, const double& z ) const
		{
			double val = 0;
			for( int n=0; n<mononomials_.size(); ++n )
			{
				val += mononomials_[n].eval(x,y,z);
			}
			return val;
		}

		double operator()( const double point[NDIM] ) const
		{
			return eval(point);
		}

		double operator()( const double& x ) const
		{
			return eval(x);
		}

		double operator()( const double& x, const double& y ) const
		{
			return eval(x,y);
		}

		double operator()( const double& x, const double& y, const double& z ) const
		{
			return eval(x,y,z);
		}

		Polynomial deriv(int dim) const
		{
			Polynomial p;
			for( int i=0; i<mononomials_.size(); ++i )
				p.add(mononomials_[i].deriv(dim));
			return p;
		}

		std::vector<Polynomial> grad() const
		{
			std::vector<Polynomial> g;
			for( int d=0; d<NDIM; ++d )
				g.push_back(deriv(d));
			return g;
		}

		static Polynomial div( const std::vector<Polynomial>& pvec )
		{
			Polynomial p;
			for( int d=0; d<NDIM; ++d )
				p += pvec[d].deriv(d);
			return p;
		}

		static Polynomial curl_z( const std::vector<Polynomial>& pvec )
		{
			return pvec[1].deriv(0) - pvec[0].deriv(1);
		}

		static std::vector<Polynomial> curl( const std::vector<Polynomial>& pvec )
		{
			std::vector<Polynomial> p(3);
			if( NDIM==2 )
			{
				p[2] = pvec[1].deriv(0) - pvec[0].deriv(1);
			}
			if( NDIM==3 )
			{
				p[0] = pvec[2].deriv(1) - pvec[1].deriv(2);
				p[1] = pvec[0].deriv(2) - pvec[2].deriv(0);
				p[2] = pvec[1].deriv(0) - pvec[0].deriv(1);
			}
			return p;
		}

		Polynomial operator+( const Polynomial& p2 )
		{
			Polynomial p = *this;
			p += p2;
			return p;
		}

		Polynomial operator-( const Polynomial& p2 )
		{
			Polynomial p = *this;
			p -= p2;
			return p;
		}

		Polynomial& operator+=( const Polynomial& other )
		{
			return addition(other,+1.);
		}

		Polynomial& operator-=( const Polynomial& other )
		{
			return addition(other,-1.);
		}

	private:
		Polynomial& addition(const Polynomial& other, double sign = 1. )
		{
			std::vector<bool> matched(other.mononomials_.size(),false);
			for( int i=0; i<mononomials_.size(); ++i )
			{
				for( int j=0; j<other.mononomials_.size(); ++j )
				{
					if( !matched[j] )
					{
						bool match=true;
						for( int d=0; d<NDIM; ++d )
						{
							if( mononomials_[i].power[d] != other.mononomials_[j].power[d] )
							{
								match = false;
								break;
							}
						}
						if( match )
						{
							matched[j] = true;
							mononomials_[i].coeff += sign * other.mononomials_[j].coeff;
						}
					}
				}
			}
			for( int j=0; j<other.mononomials_.size(); ++j )
			{
				if( !matched[j] )
				{
					add(other.mononomials_[j]);
				}
			}
			return *this;
		}



	private:
		std::vector< monomial_type > mononomials_;
	};

	class Element : public eckit::Owned
	{
	public:
		typedef eckit::SharedPtr<Element      > Ptr;
		typedef std::vector< Element::Ptr > Vector;

	public:

		int N() const { return N_; }

		const double* nodes() const { return nodes_.data(); }

	protected:

		int N_;

		std::vector<double> nodes_;

	};

	class Point : public Element
	{
	public:
		Point()
		{
			N_ = 1;
		}
	};


	class QuadP1 : public Element
	{
	public:
		QuadP1()
		{
			N_ = 4;
			int nodes_init[] = {
				-1., -1.,
				 1., -1.,
				 1.,  1.,
				-1.,  1. };
			nodes_.assign(nodes_init, nodes_init+N_);
		}
	};

	class TriagP1 : public Element
	{
	public:
		TriagP1()
		{
			N_ = 3;
			int nodes_init[] = {
				 0., 0.,
				 1., 0.,
				 0., 1. };
			nodes_.assign(nodes_init, nodes_init+N_);
		}
	};

	class LineP0 : public Element
	{
	public:
		LineP0()
		{
			N_ = 1;
			int nodes_init[] = { 0. };
			nodes_.assign(nodes_init, nodes_init+N_);
		}
	};

	class LineP1 : public Element
	{
	public:
		LineP1()
		{
			N_ = 2;
			int nodes_init[] = { -1., 1. };
			nodes_.assign(nodes_init, nodes_init+N_);
		}
	};

	class LineP2 : public Element
	{
	public:
		LineP2()
		{
			N_ = 3;
			int nodes_init[] = { -1., 1., 0. };
			nodes_.assign(nodes_init, nodes_init+N_);
		}
	};

	class Structured1D: public Element
	{
	public:
		Structured1D(int N)
		{
			N_ = N;
		}
	};


	class Nodes : public eckit::Owned
	{
	public:
		typedef eckit::SharedPtr<Nodes> Ptr;
	public:
		Nodes(int ngp, int nlev)
		{
			ngp_  = ngp;
			nlev_ = nlev;
		}
		std::vector< double > field;
		int ngp_;
		int nlev_;
	};

	class Elements
	{
	public:
//		Elements(int nlev)
//		{
//			nlev_ = nlev;
//		}
//		Elements& allocate_dof(int nnodes)
//		{
//			field_.resize(nnodes*nlev_);
//			return *this;
//		}
//		Elements& allocate_dof()
//		{
//			int dof=0;
//			for( int etype=0; etype<elements_.size(); ++etype )
//				dof += nelem_[etype]*elements_[etype]->N();
//			dof *= nlev_;
//			field_.resize(dof);
//			return *this;
//		}
		Elements& add(const std::string& element_type, int nelem)
		{
			if( element_type == "QuadP1" )   elements_.push_back( Element::Ptr( new QuadP1 ) );
			if( element_type == "TriagP1" )  elements_.push_back( Element::Ptr( new TriagP1 ) );
			if( element_type == "LineP0" )   elements_.push_back( Element::Ptr( new LineP0 ) );
			if( element_type == "LineP1" )   elements_.push_back( Element::Ptr( new LineP1 ) );
			if( element_type == "LineP2" )   elements_.push_back( Element::Ptr( new LineP2 ) );
			index_[element_type] = nelem_.size();
			nelem_.push_back(nelem);
			return *this;
		}
		Elements& nodes( const Nodes::Ptr& nodes )
		{
			nodes_ = nodes;
			return *this;
		}

		std::map<std::string,int> index_;
		Element::Vector elements_;
		Nodes::Ptr nodes_;
		std::vector<int> nelem_;
		std::vector< atlas::Array<double> > P0_fields;
	};

	class Column // Really is function space
	{
	public:
		Column(const Element::Ptr& horizontal, const Element::Ptr& vertical):
			horizontal_(horizontal),
			vertical_(vertical)
		{
		}
		int N() const { return ngp()*nlev(); }
		int ngp() const { return horizontal_->N(); }
		int nlev() const { return vertical_->N(); }
	protected:
		Element::Ptr horizontal_;
		Element::Ptr vertical_;
	};




	template< int NDIM >
	std::vector< Polynomial<NDIM> > polynomial_basis(int order, double* points, int npts)
	{
		std::vector< Polynomial<NDIM> > basis(npts);

		Matrix<double> vandermonde(npts,npts);
		Matrix<double> coefficients(npts,npts);
		Matrix<int>    powers(npts,NDIM);
		Matrix<double>::Proxy pts(points,npts,NDIM);
		if (NDIM==1)
		{
			size_t n(0);
      for(size_t k=0; k<=order; ++k)
      {
        powers(n,0) = k;
        ++n;
      }
      for (n=0; n<npts; ++n)
      {
        double x = pts(n,0);
        for (size_t p=0; p<npts; ++p)
        {
          double x_p = ( powers(p,0) == 0. ? 1.: std::pow(x,powers(p,0)) );
          vandermonde(n,p) = x_p;
        }
      }
    }
		else if (NDIM==2)
		{
			size_t n(0);
			for(size_t l=0; l<=order; ++l) {
				for(size_t k=0; k<=order; ++k) {
					powers(n,0) = k;
					powers(n,1) = l;
					++n;
				}
			}
			for (n=0; n<npts; ++n)
			{
				double x = pts(n,0);
				double y = pts(n,1);
				for (size_t p=0; p<npts; ++p)
				{
					double x_p = ( powers(p,0) == 0. ? 1.: std::pow(x,powers(p,0)) );
					double y_p = ( powers(p,1) == 0. ? 1.: std::pow(y,powers(p,1)) );
					vandermonde(n,p) = x_p * y_p;
				}
			}
		}
		std::cout << "pts = \n" << pts << std::endl;
		std::cout << "powers = \n" << powers << std::endl;
		std::cout << "vandermonde = \n" << vandermonde << std::endl;
		coefficients = vandermonde.inverse().transpose();
		std::cout << "coefficients = \n" << coefficients << std::endl;
		int pwr[3];
		for( size_t n=0; n<npts; ++n )
		{
			for( size_t i=0; i<npts; ++i )
			{
				for( int d=0; d<NDIM; ++d )
					pwr[d] = powers(i,d);
				basis[n].add( Mononomial<NDIM>( coefficients(n,i), pwr ) );
			}
		}
		return basis;
	}

} // namespace test
} // namespace atlas

using namespace atlas::test;

BOOST_AUTO_TEST_CASE( test_functionspace )
{

	atlas::mpi::init();

	Element::Ptr point( new Point );
	Element::Ptr quad( new QuadP1 );
	Element::Ptr levels( new Structured1D(100) );

	Column quad_column(quad,levels);
	Column point_column(point,levels);

	//std::cout << quad->N() << std::endl;
	//std::cout << levels->N() << std::endl;
	DEBUG_VAR(quad_column.N());
	DEBUG_VAR(quad_column.nlev());
	DEBUG_VAR(quad_column.ngp());

	DEBUG_VAR(point_column.N());
	DEBUG_VAR(point_column.nlev());
	DEBUG_VAR(point_column.ngp());

	Nodes::Ptr nodes( new Nodes(7,10) );

	Elements elements;
	elements.nodes( nodes );
	elements.add("QuadP1",2).add("TriagP1",1);

	Elements edges;
	edges.nodes( nodes );
	edges.add("LineP0",9);

	typedef Polynomial<2> polynomial_type;
	typedef polynomial_type::monomial_type monomial_type;


	polynomial_type p;

	p.add( monomial_type( 3., 0,0 ) );
	p.add( monomial_type( 1., 1,0 ) );
	p.add( monomial_type( 1., 2,0 ) );
	DEBUG_VAR(p(2,2));

	polynomial_type dpdx = p.deriv(0);
	DEBUG_VAR(dpdx(2,2));

	polynomial_type dpdy = p.deriv(1);
	DEBUG_VAR(dpdy(2,2));

	std::vector<polynomial_type> grad = p.grad();

	std::vector<polynomial_type> pvec(2,p);
	polynomial_type div = polynomial_type::div(pvec);
	DEBUG_VAR(div(2,2));

	DEBUG_VAR( polynomial_type::curl_z(grad)(2,2) );

	typedef std::vector< Polynomial<1> > PolynomialBasis1D;
	typedef std::vector< Polynomial<2> > PolynomialBasis2D;


	Matrix<int> m;
//	BOOST_CHECK( m.data_ == NULL );
//	BOOST_CHECK( m.nr_ == 0 );
//	BOOST_CHECK( m.nc_ == 0 );
	m.resize(2,3);

	BOOST_CHECK_EQUAL(m.size(),6);
	BOOST_CHECK_EQUAL(m.rows(),2);
	BOOST_CHECK_EQUAL(m.cols(),3);

	m(0,0) = 0;   m(0,1) = 2;   m(0,2) = 4;
	m(1,0) = 1;   m(1,1) = 3;   m(1,2) = 5;

	m *= 10;

	Matrix<int> b = m + m;

	std::cout << "m = \n" << b << std::endl;

	double line_pts[] = {-1. , 1.};
	int line_npts = 2;
	PolynomialBasis1D line_basis = polynomial_basis<1>(1,line_pts,line_npts);
	for( int n=0; n<line_npts; ++n)
		DEBUG( n << ":  " << line_basis[n](-1.) << "  " << line_basis[n](+1.) );

	double quad_pts[] = {-1,   1.,   1.,  -1,
											 -1., -1.,   1.,   1.};
	int quad_npts = 4;
	PolynomialBasis2D quad_basis = polynomial_basis<2>(1,quad_pts,quad_npts);

	for( int n=0; n<quad_npts; ++n)
		DEBUG( n << ":  " << quad_basis[n](-1.,-1.) << "  " << quad_basis[n](1.,-1.) << " " << quad_basis[n](1.,1.) << " " << quad_basis[n](-1.,1.) );


	RowVector<int> vr(3);
	vr(0) = 0; vr(1) = 1; vr(2) = 2;

	ColVector<int> vc(3);
	vc(0) = 0; vc(1) = 1; vc(2) = 2;

	std::cout << vr << "\n" << vc << "\n\n" << vr*vc << std::endl;



	mpi::finalize();
}
