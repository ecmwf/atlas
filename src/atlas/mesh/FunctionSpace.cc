/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#include <cassert>
#include <iostream>
#include <sstream>
#include <limits>
#include <eckit/exception/Exceptions.h>
#include "atlas/atlas_defines.h"
#include "atlas/mesh/FunctionSpace.h"
#include "atlas/mesh/Field.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/util/Debug.h"

#ifdef HAVE_FORTRAN_NUMBERING
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif


namespace atlas {

FunctionSpace::FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<int>& shape ) :
	name_(name),
	shape_(shape),
	gather_scatter_(new mpl::GatherScatter()),
	halo_exchange_(new mpl::HaloExchange()),
	checksum_(new mpl::Checksum()),
	mesh_(NULL)
{
	//std::cout << "C++ : shape Constructor" << std::endl;
	dof_ = 1;
	size_t extsize = shape_.size();
	shapef_.resize(extsize);
	for (size_t i=0; i<extsize; ++i)
	{
		shapef_[extsize-1-i] = shape_[i];
		if( shape_[i] != Field::UNDEF_VARS )
			dof_ *= shape_[i];
	}
	glb_dof_ = dof_;
}

FunctionSpace::~FunctionSpace()
{
//	std::cout << "shape Destructor ("<<name_<<")" << std::endl;
	index_.clear();
	for( size_t f=0; f<fields_.size(); ++f )
		if( fields_[f] ) delete(fields_[f]);
	fields_.clear();
}

void FunctionSpace::resize(const std::vector<int>& shape)
{
	if (shape.size() != shape_.size() )
		throw eckit::BadParameter("Cannot resize shape: shape sizes don't match.",Here());

	size_t extsize = shape_.size();

	for (size_t i=1; i<extsize; ++i)
	{
		if (shape[i] != shape_[i])
			throw eckit::BadParameter("Only the first extent can be resized for now!",Here());
	}

	shape_ = shape;
	shapef_.resize(extsize);
	for (size_t i=0; i<extsize; ++i)
	{
		shapef_[extsize-1-i] = shape_[i];
	}

	dof_ = 1;
	for (size_t i=0; i<extsize; ++i)
	{
		if( shape_[i] != Field::UNDEF_VARS )
			dof_ *= shape_[i];
	}

	for( int f=0; f<fields_.size(); ++f)
	{
		std::vector< int > field_shape(extsize);
		for (size_t i=0; i<extsize; ++i)
		{
			if( shape_[i] == Field::UNDEF_VARS )
				field_shape[i] = fields_[f]->nb_vars();
			else
				field_shape[i] = shape_[i];
		}
		fields_[f]->allocate(field_shape);
	}
}

template <>
FieldT<double>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
	if( has_field(name) )
	{
		std::ostringstream msg; msg << "field with name " << name << "already exists" << std::endl;
		throw eckit::Exception( msg.str(), Here() );
	}

	index_[name] = fields_.size();
	FieldT<double>* field = new FieldT<double>(name,nb_vars,*this);
	fields_.push_back( field );

	size_t rank = shape_.size();
	std::vector< int > field_shape(rank);
	for (size_t i=0; i<rank; ++i)
	{
		if( shape_[i] == Field::UNDEF_VARS )
			field_shape[i] = field->nb_vars();
		else
			field_shape[i] = shape_[i];
	}

	field->allocate(field_shape);
	return *field;
}

template <>
FieldT<float>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
	if( has_field(name) )
	{
		std::ostringstream msg; msg << "field with name " << name << "already exists" << std::endl;
		throw eckit::Exception( msg.str(), Here() );
	}

	index_[name] = fields_.size();
	FieldT<float>* field = new FieldT<float>(name,nb_vars,*this);
	fields_.push_back( field );

	size_t rank = shape_.size();
	std::vector< int > field_shape(rank);
	for (size_t i=0; i<rank; ++i)
	{
		if( shape_[i] == Field::UNDEF_VARS )
			field_shape[i] = field->nb_vars();
		else
			field_shape[i] = shape_[i];
	}

	field->allocate(field_shape);
	return *field;
}

template <>
FieldT<int>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
	if( has_field(name) )
	{
		std::ostringstream msg; msg << "field with name " << name << "already exists" << std::endl;
		throw eckit::Exception( msg.str(), Here() );
	}

	index_[name] = fields_.size();
	FieldT<int>* field = new FieldT<int>(name,nb_vars,*this);
	fields_.push_back( field );


	size_t rank = shape_.size();
	std::vector< int > field_shape(rank);
	for (size_t i=0; i<rank; ++i)
	{
		if( shape_[i] == Field::UNDEF_VARS )
			field_shape[i] = field->nb_vars();
		else
			field_shape[i] = shape_[i];
	}

	field->allocate(field_shape);
	return *field;
}

void FunctionSpace::remove_field(const std::string& name)
{
	if( has_field(name) )
	{
		delete( fields_[ index_.at(name) ] );
		fields_[ index_.at(name) ] = 0;
		index_.erase(name);
	}
	else
	{
		std::stringstream msg;
		msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
		throw eckit::OutOfRange(msg.str(),Here());
	}
}

Field& FunctionSpace::field( size_t idx ) const
{
	assert( idx < fields_.size() );
	return *fields_[ idx ];
}

Field& FunctionSpace::field(const std::string& name) const
{
	if( has_field(name) )
	{
		return *fields_[ index_.at(name) ];
	}
	else
	{
		std::stringstream msg;
		msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
		throw eckit::OutOfRange(msg.str(),Here());
	}
}

template<>
	FieldT<double> &FunctionSpace::field(const std::string& name) const
{
	if( has_field(name) )
	{
		return *dynamic_cast< FieldT<double>* >(fields_[ index_.at(name) ]);
	}
	else
	{
		std::stringstream msg;
		msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
		throw eckit::OutOfRange(msg.str(),Here());
	}
}

template<>
	FieldT<float> &FunctionSpace::field(const std::string& name) const
{
	if( has_field(name) )
	{
		return *dynamic_cast< FieldT<float>* >(fields_[ index_.at(name) ]);
	}
	else
	{
		std::stringstream msg;
		msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
		throw eckit::OutOfRange(msg.str(),Here());
	}
}

template<>
	FieldT<int> &FunctionSpace::field(const std::string& name) const
{
	if( has_field(name) )
	{
		return *dynamic_cast< FieldT<int>* >(fields_[ index_.at(name) ]);
	}
	else
	{
		std::stringstream msg;
		msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
		throw eckit::OutOfRange(msg.str(),Here());
	}
}

void FunctionSpace::parallelise(const int part[], const int remote_idx[], const int glb_idx[], int parsize)
{
	halo_exchange_->setup(part,remote_idx,REMOTE_IDX_BASE,parsize);
	gather_scatter_->setup(part,remote_idx,REMOTE_IDX_BASE,glb_idx,-1,parsize);
	checksum_->setup(part,remote_idx,REMOTE_IDX_BASE,glb_idx,-1,parsize);
	glb_dof_ = gather_scatter_->glb_dof();
	for( int b=shapef_.size()-2; b>=0; --b)
	{
		if( shapef_[b] != Field::UNDEF_VARS )
			glb_dof_ *= shapef_[b];
	}
}

void FunctionSpace::parallelise(FunctionSpace& other_shape)
{
	halo_exchange_ = mpl::HaloExchange::Ptr( &other_shape.halo_exchange() );
	gather_scatter_ = mpl::GatherScatter::Ptr( &other_shape.gather_scatter() );
}

void FunctionSpace::parallelise()
{
	if( name() == "nodes" || name() == "edges" )
	{
		FieldT<int>& ridx = field<int>("remote_idx");
		FieldT<int>& part = field<int>("partition");
		FieldT<int>& gidx = field<int>("glb_idx");
		parallelise(part.data(),ridx.data(),gidx.data(),part.size());
	}
	else
	{
		NOTIMP;
	}
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FunctionSpace* atlas__FunctionSpace__new (char* name, char* shape_func, int shape[], int shape_size) {
	std::vector<int> shape_vec(shape,shape+shape_size);
	return new FunctionSpace( std::string(name), std::string(shape_func), shape_vec );
}

int atlas__FunctionSpace__dof (FunctionSpace* This) {
	return This->dof();
}

int atlas__FunctionSpace__glb_dof (FunctionSpace* This) {
	return This->glb_dof();
}

void atlas__FunctionSpace__create_field_double (FunctionSpace* This, char* name, int nb_vars) {
	This->create_field<double>( std::string(name), nb_vars );
}

void atlas__FunctionSpace__create_field_float (FunctionSpace* This, char* name, int nb_vars) {
	This->create_field<float>( std::string(name), nb_vars );
}

void atlas__FunctionSpace__create_field_int (FunctionSpace* This, char* name, int nb_vars) {
	This->create_field<int>( std::string(name), nb_vars );
}

void atlas__FunctionSpace__remove_field (FunctionSpace* This, char* name ) {
	This->remove_field( std::string(name) );
}

int atlas__FunctionSpace__has_field (FunctionSpace* This, char* name) {
	return This->has_field( std::string(name) );
}

const char* atlas__FunctionSpace__name (FunctionSpace* This) {
	return This->name().c_str();
}

void atlas__FunctionSpace__shapef (FunctionSpace* This, int* &shape, int &rank) {
	shape = const_cast<int*>(&(This->shapef()[0]));
	rank = This->shapef().size();
}

Field* atlas__FunctionSpace__field (FunctionSpace* This, char* name) {
	return &This->field( std::string(name) );
}

void atlas__FunctionSpace__parallelise (FunctionSpace* This) {
	This->parallelise();
}

void atlas__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size) {
	This->halo_exchange(field_data,field_size);
}

void atlas__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size) {
	This->halo_exchange(field_data,field_size);
}

void atlas__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size) {
	This->halo_exchange(field_data,field_size);
}

void atlas__FunctionSpace__gather_int (FunctionSpace* This, int field_data[], int field_size, int glbfield_data[], int glbfield_size) {
	This->gather(field_data,field_size, glbfield_data,glbfield_size);
}

void atlas__FunctionSpace__gather_float (FunctionSpace* This, float field_data[], int field_size, float glbfield_data[], int glbfield_size) {
	This->gather(field_data,field_size, glbfield_data,glbfield_size);
}

void atlas__FunctionSpace__gather_double (FunctionSpace* This, double field_data[], int field_size, double glbfield_data[], int glbfield_size) {
	This->gather(field_data,field_size, glbfield_data,glbfield_size);
}

mpl::HaloExchange* atlas__FunctionSpace__halo_exchange (FunctionSpace* This) {
	return &This->halo_exchange();
}

mpl::GatherScatter* atlas__FunctionSpace__gather (FunctionSpace* This) {
	return &This->gather_scatter();
}

mpl::Checksum* atlas__FunctionSpace__checksum (FunctionSpace* This) {
	return &This->checksum();
}

void atlas__FunctionSpace__delete (FunctionSpace* This) {
	delete This;
}
// ------------------------------------------------------------------

} // namespace atlas

