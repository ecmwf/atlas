#include "atlas_trace.h"

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/CodeLocation.h"

namespace atlas {

//static std::vector<std::string> Labels( int size, const char* labels[] ) {
//    std::vector<std::string> _labels;
//    _labels.reserve( size );
//    for ( int i = 0; i < size; ++i ) {
//        _labels.emplace_back( labels[i] );
//    }
//    return _labels;
//}

extern "C" {

Trace* new_atlas_Trace( const char* file, int line, const char* title ) {
    return new Trace( CodeLocation( file, line, nullptr, true ), std::string( title ) );
}

Trace* new_atlas_Trace_labels_1( const char* file, int line, const char* title, const char* label1 ) {
    std::vector<std::string> labels{label1};
    return new Trace( CodeLocation( file, line, nullptr, true ), std::string( title ), labels );
}
Trace* new_atlas_Trace_labels_2( const char* file, int line, const char* title, const char* label1,
                                 const char* label2 ) {
    std::vector<std::string> labels{label1, label2};
    return new Trace( CodeLocation( file, line, nullptr, true ), std::string( title ), labels );
}
Trace* new_atlas_Trace_labels_3( const char* file, int line, const char* title, const char* label1, const char* label2,
                                 const char* label3 ) {
    std::vector<std::string> labels{label1, label2, label3};
    return new Trace( CodeLocation( file, line, nullptr, true ), std::string( title ), labels );
}
Trace* new_atlas_Trace_labels_4( const char* file, int line, const char* title, const char* label1, const char* label2,
                                 const char* label3, const char* label4 ) {
    std::vector<std::string> labels{label1, label2, label3, label4};
    return new Trace( CodeLocation( file, line, nullptr, true ), std::string( title ), labels );
}
Trace* new_atlas_Trace_labels_5( const char* file, int line, const char* title, const char* label1, const char* label2,
                                 const char* label3, const char* label4, const char* label5 ) {
    std::vector<std::string> labels{label1, label2, label3, label4, label5};
    return new Trace( CodeLocation( file, line, nullptr, true ), std::string( title ), labels );
}

void delete_atlas_Trace( Trace* This ) {
    ASSERT( This != nullptr );
    delete This;
    This = nullptr;
}

void atlas_Trace__start( Trace* This ) {
    ASSERT( This != nullptr );
    This->start();
}

void atlas_Trace__stop( Trace* This ) {
    ASSERT( This != nullptr );
    This->stop();
}

void atlas_Trace__pause( Trace* This ) {
    ASSERT( This != nullptr );
    This->pause();
}

void atlas_Trace__resume( Trace* This ) {
    ASSERT( This != nullptr );
    This->resume();
}

int atlas_Trace__running( Trace* This ) {
    ASSERT( This != nullptr );
    return This->running();
}

double atlas_Trace__elapsed( Trace* This ) {
    ASSERT( This != nullptr );
    return This->elapsed();
}

}  // extern C

}  // namespace atlas
