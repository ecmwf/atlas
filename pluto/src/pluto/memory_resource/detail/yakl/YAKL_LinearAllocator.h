
#pragma once
// Included by YAKL_Gator.h

__YAKL_NAMESPACE_WRAPPER_BEGIN__
namespace yakl {


  // This class encapsulates a single "pool"
  /** @private */
  class LinearAllocator {
  public:

    // Describes a single allocation entry
    struct AllocNode {
      size_t start;        // Offset of this allocation in "blocks"
      size_t length;       // Length of this allocation in "blocks"
      char const * label;  // Label for this allocation
      AllocNode() {
        this->start  = 0;
        this->length = 0;
        this->label  = "";
      }
      AllocNode( size_t start , size_t length , char const * label ) {
        this->start  = start;
        this->length = length;
        this->label  = label;
      }
    };

    std::string                           pool_name;
    void                                  *pool;     // Raw pool pointer
    size_t                                nBlocks;   // Number of blocks in the pool
    unsigned                              blockSize; // Size of each block in bytes
    unsigned                              blockInc;  // Number of size_t in each block
    std::vector<AllocNode>                allocs;    // List of allocations
    std::function<void *( size_t )>       mymalloc;  // allocation function
    std::function<void( void * )>         myfree;    // free function
    std::function<void( void *, size_t )> myzero;    // zero function


    LinearAllocator() { nullify(); }


    LinearAllocator( size_t                                bytes ,
                     unsigned                              blockSize = 16*sizeof(size_t) ,
                     std::function<void *( size_t )>       mymalloc  = [] (size_t bytes) -> void * { return ::malloc(bytes); } ,
                     std::function<void( void * )>         myfree    = [] (void *ptr) { ::free(ptr); } ,
                     std::function<void( void *, size_t )> myzero    = [] (void *ptr, size_t bytes) {} ,
                     std::string                           pool_name = "Gator" ,
                     std::string                           error_message_out_of_memory = "" ) {
      nullify();

      #ifdef YAKL_VERBOSE
        if (bytes >= 1024*1024*1024) {
          verbose_inform(std::string("Creating pool of ")+std::to_string(bytes/1024./1024./1024.)+" GB" , pool_name);
        } else if (bytes >= 1024*1024) {
          verbose_inform(std::string("Creating pool of ")+std::to_string(bytes/1024./1024.      )+" MB" , pool_name);
        } else if (bytes >= 1024) {
          verbose_inform(std::string("Creating pool of ")+std::to_string(bytes/1024.            )+" KB" , pool_name);
        } else {
          verbose_inform(std::string("Creating pool of ")+std::to_string(bytes                  )+" B"  , pool_name);
        }
      #endif
      if (blockSize%(2*sizeof(size_t)) != 0) {
        std::cerr << "ERROR: Pool labeled \"" << pool_name << "\" -> LinearAllocator:" << std::endl;
        die("Error: LinearAllocator blockSize must be a multiple of 2*sizeof(size_t)");
      }
      this->blockSize = blockSize;
      this->blockInc  = blockSize / sizeof(size_t); // In two routines, I cast void * to size_t * for arithmetic
                                                    // Therefore, blockInc is the number of size_t in a block
      this->nBlocks   = (bytes-1) / blockSize + 1;
      this->mymalloc  = mymalloc;
      this->myfree    = myfree  ;
      this->myzero    = myzero  ;
      this->pool      = mymalloc( poolSize() );
      this->allocs    = std::vector<AllocNode>();
      this->allocs.reserve(128);  // Make sure there is initial room for 128 entries
      this->pool_name = pool_name;
      if (pool == nullptr) {
        std::cerr << "ERROR: Pool labeled \"" << pool_name << "\" -> LinearAllocator:" << std::endl;
        std::cerr << "Could not create pool of size " << bytes << " bytes (" << bytes/1024./1024./1024. << " GB)."
                  << "\nYou have run out of memory." << std::endl;
        std::cerr << "When individual variables consume sizable percentages of a pool, memory gets segmented, and "
                  << "the pool space isn't used efficiently. \nLarger pools will improve that. So try increasing the "
                  << "size of the initial pool and maybe the grow size as well. \nIn the extreme, you could create "
                  << "an initial pool that consumes most of the avialable memory. \nIf that still doesn't work, then "
                  << "it sounds like you're choosing a problem size that's too large for the number of compute "
                  << "nodes you're using.\n";
        die( error_message_out_of_memory );
      }
      this->myzero( pool , poolSize() );
    }


    // Allow the pool to be moved, but not copied
    LinearAllocator( LinearAllocator && rhs) {
      this->pool      = rhs.pool     ;
      this->nBlocks   = rhs.nBlocks  ;
      this->blockSize = rhs.blockSize;
      this->blockInc  = rhs.blockInc ;
      this->allocs    = rhs.allocs   ;
      this->mymalloc  = rhs.mymalloc ;
      this->myfree    = rhs.myfree   ;
      this->myzero    = rhs.myzero   ;
      this->pool_name = rhs.pool_name;
      rhs.nullify();
    }


    LinearAllocator &operator =( LinearAllocator && rhs) {
      if (this == &rhs) { return *this; }
      this->finalize();
      this->pool      = rhs.pool     ;
      this->nBlocks   = rhs.nBlocks  ;
      this->blockSize = rhs.blockSize;
      this->blockInc  = rhs.blockInc ;
      this->allocs    = rhs.allocs   ;
      this->mymalloc  = rhs.mymalloc ;
      this->myfree    = rhs.myfree   ;
      this->myzero    = rhs.myzero   ;
      this->pool_name = rhs.pool_name;
      rhs.nullify();
      return *this;
    }


    LinearAllocator( LinearAllocator const &rhs ) = delete;


    LinearAllocator &operator=( LinearAllocator const &rhs ) = delete;


    ~LinearAllocator() {
      if (pool != nullptr) {
        verbose_inform("Destroying pool" , pool_name);
      }
      finalize();
    }


    void nullify() {
      this->pool      = nullptr;
      this->nBlocks   = 0;
      this->blockSize = 0;
      this->blockInc  = 0;
      this->allocs    = std::vector<AllocNode>();
      this->mymalloc  = [] (size_t bytes) -> void * { return ::malloc(bytes); };
      this->myfree    = [] (void *ptr) { ::free(ptr); };
      this->myzero    = [] (void *ptr, size_t bytes) {};
    }


    void finalize() {
      if (allocs.size() != 0) {
        #if defined(YAKL_DEBUG)
          std::cerr << "WARNING: Pool labeled \"" << pool_name << "\" -> LinearAllocator:" << std::endl;
          std::cerr << "WARNING: Not all allocations were deallocated before destroying this pool.\n" << std::endl;
          printAllocsLeft();
          std::cerr << "This probably won't end well, but carry on.\n" << std::endl;
        #endif
      }
      if (this->pool != nullptr) { myfree( this->pool ); }
      nullify();
    }


    // Mostly for debug purposes. Print all existing allocations
    void printAllocsLeft() {
      if (allocs.size() != 0) {
        std::cerr << "The following allocations have not been deallocated:" << std::endl;
        for (int i=0; i < allocs.size(); i++) {
          std::cerr << "*** Label: " << allocs[i].label
                    << "  ;  size: " << allocs[i].length*blockSize
                    << " bytes  ;  offset: " << allocs[i].start*blockSize
                    << " bytes  ;  ptr: " << getPtr(allocs[i].start) << std::endl;
        }
      }
    }


    // Allocate the requested number of bytes if there is room for it.
    // If there isn't room or bytes == 0, then nullptr is returned.
    // Otherwise, the correct pointer is returned
    void * allocate(size_t bytes, char const * label="") {
      if (bytes == 0) {
        return nullptr;
      }
      size_t blocksReq = (bytes-1)/blockSize + 1; // Number of blocks needed for this allocation
      // If there are no allocations, then place this allocaiton at the beginning.
      if (allocs.empty()) {
        if (nBlocks >= blocksReq) {
          allocs.push_back( AllocNode( (size_t) 0 , blocksReq , label ) );
          return pool;
        }
      } else {
        // Look for room before the first allocation
        if ( allocs.front().start >= blocksReq ) {
          allocs.insert( allocs.begin() , AllocNode( 0 , blocksReq , label ) );
          return getPtr(allocs[0].start);
        }

        // Loop through the allocations and look for free space between this and the next
        for (int i=0; i < allocs.size()-1; i++) {
          if ( allocs[i+1].start - (allocs[i].start + allocs[i].length) >= blocksReq ) {
            allocs.insert( allocs.begin()+i+1 , AllocNode( allocs[i].start+allocs[i].length , blocksReq , label ) );
            return getPtr(allocs[i+1].start);
          }
        }

        // Look for room after the last allocation
        if ( nBlocks - (allocs.back().start + allocs.back().length) >= blocksReq ) {
          allocs.push_back( AllocNode( allocs.back().start + allocs.back().length , blocksReq , label ) );
          return getPtr(allocs.back().start);
        }
      }

      // Return nullptr if there was no room for the allocation
      // If the caller used "iGotRoom", then this should never actually happen
      return nullptr;
    };


    // Free the requested pointer
    // Returns the number of bytes in the allocation being freed
    size_t free(void *ptr, char const * label = "") {
      for (int i=allocs.size()-1; i >= 0; i--) {
        if (ptr == getPtr(allocs[i].start)) {
          size_t bytes = allocs[i].length*blockSize;
          allocs.erase(allocs.begin()+i);
          return bytes;
        }
      }
      std::cerr << "ERROR: Pool labeled \"" << pool_name << "\" -> LinearAllocator:" << std::endl;
      std::cerr << "Trying to free an invalid pointer.\n";
      die("This means you have either already freed the pointer, or its address has been corrupted somehow.");
      return 0;
    };


    // Determine if there is room for an allocation of the requested number of bytes
    bool iGotRoom( size_t bytes ) const {
      size_t blocksReq = (bytes-1)/blockSize + 1; // Number of blocks needed for this allocation

      if (allocs.empty()) {
        if (nBlocks >= blocksReq) { return true; }
      } else {
        // Look for room before the first allocation
        if ( allocs.front().start >= blocksReq ) { return true; }

        // Loop through the allocations and look for free space between this and the next
        for (int i=0; i < allocs.size()-1; i++) {
          if ( allocs[i+1].start - (allocs[i].start + allocs[i].length) >= blocksReq ) { return true; }
        }

        // Look for room after the last allocation
        if ( nBlocks - (allocs.back().start + allocs.back().length) >= blocksReq ) { return true; }
      }

      return false;
    }


    // Determine if the requested pointer belongs to this pool
    bool thisIsMyPointer(void *ptr) const {
      long long offset = ( (size_t *) ptr - (size_t *) pool ) / blockInc;
      return (offset >= 0 && offset <= nBlocks-1);
    }


    bool initialized() const { return pool != nullptr; }


    size_t poolSize() const { return nBlocks*blockSize; }


    size_t numAllocs() const { return allocs.size(); }


    // Transform a block index into a memory pointer
    void * getPtr( size_t blockIndex ) const {
      return (void *) ( ( (size_t *) pool ) + blockIndex*blockInc );
    }


    void die(std::string str="") {
      std::cerr << str << std::endl;
      throw std::runtime_error(str);
    }


  };

}
__YAKL_NAMESPACE_WRAPPER_END__


