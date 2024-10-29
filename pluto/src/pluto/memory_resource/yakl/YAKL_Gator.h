/**
 * @file
 * The YAKL pool allocator class
 */

#pragma once
// Included by YAKL.h

#include "YAKL_LinearAllocator.h"

__YAKL_NAMESPACE_WRAPPER_BEGIN__
namespace yakl {

  /** @brief YAKL Pool allocator class.
    *
    * Growable pool allocator for efficient frequent allocations and deallocations.
    * User determines allocation, free, initial pool size, additional pool size, and other pool allocator
    * characteristics upon initiailization.
    *
    * Once existing pools run out of memory, additional pools are create. Each pool is based on a simple
    * linear search for free slots that is as efficient in memory usage as possible. While search time is
    * linear rather than log in complexity, allocations and frees are typically overlapped with kernel execution.
    *
    * Gator objects **are thread safe** for allocate() and free() calls.
    */
  class Gator {
  protected:
    /** @private */
    std::string                           pool_name;
    /** @private */
    std::list<LinearAllocator>            pools;    // The pools managed by this class
    /** @private */
    std::function<void *( size_t )>       mymalloc; // allocation function
    /** @private */
    std::function<void( void * )>         myfree;   // free function
    /** @private */
    std::function<void( void *, size_t )> myzero;   // zero function
    /** @private */
    size_t growSize;   // Amount by which the pool grows in bytes
    /** @private */
    size_t blockSize;  // Minimum allocation size
    /** @private */
    std::string error_message_cannot_grow;
    /** @private */
    std::string error_message_out_of_memory;
    /** @private */
    size_t high_water_mark;  // memory high water mark
    /** @private */
    size_t bytes_currently_allocated;  // memory high water mark

    /** @private */
    std::mutex mtx1;    // Internal mutex used to protect alloc and free in threaded regions
    std::mutex mtx2;    // Internal mutex used to protect alloc and free in threaded regions

    /** @private */
    void die(std::string str="") {
      std::cerr << str << std::endl;
      throw std::runtime_error(str);
    }


    struct WaitEntry {
      std::string label;
      void *ptr;
      std::vector<Event> events;
    };

    std::vector<WaitEntry> waiting_entries;
    std::vector<Event>     waiting_events;


  public:

    /** @brief Please use the init() function to specify parameters, not the constructor. */
    Gator() {
    }


    // No moves allowed
    /** @brief A Gator object may be moved but not copied */
    Gator            (      Gator && );
    /** @brief A Gator object may be moved but not copied */
    Gator &operator= (      Gator && );
    /** @brief A Gator object may be moved but not copied */
    Gator            (const Gator &  ) = delete;
    /** @brief A Gator object may be moved but not copied */
    Gator &operator= (const Gator &  ) = delete;


    /** @brief All pools are automatically finalized when a Gator object is destroyed. */
    ~Gator() {
      finalize();
    }


    /** @brief Initialize the pool.
      *
      * @param myalloc The allocator to use when creating the initial and additional pools
      * @param myfree  The deallocator to use when destroying all pools upon the Gator object destruction.
      * @param myzero  [NOT CURRENTLY USED] The function to use to memset the memory to a value
      *                (presumably zero) after initial allocationa and an entry free
      * @param initialSize Size of the initial pool in bytes
      * @param growSize    Size of additional pools once the previous pools run out of memory
      * @param blockSize   Size of an individual block (the smallest possible entry size. The entry
      *                    size must increment by this size as well.
      * @param pool_name   The label for this pool used in debugging
      * @param error_message_out_of_memory The string printed when a pool allocation fails
      * @param error_message_cannot_grow   The string printed when an allocation is requested that overflows
      *                                    initial memory and exceeds the size of additional pools
      */
    void init(std::function<void *( size_t )>       mymalloc  = [] (size_t bytes) -> void * { return ::malloc(bytes); },
              std::function<void( void * )>         myfree    = [] (void *ptr) { ::free(ptr); }                        ,
              std::function<void( void *, size_t )> myzero    = [] (void *ptr, size_t bytes) {}                        ,
              size_t initialSize                              = 1024*1024*1024                                         ,
              size_t growSize                                 = 1024*1024*1024                                         ,
              size_t blockSize                                = 16*sizeof(size_t)                                      ,
              std::string pool_name                           = "Gator"                                                ,
              std::string error_message_out_of_memory         = ""                                                     ,
              std::string error_message_cannot_grow           = "" ) {
      this->mymalloc  = mymalloc ;
      this->myfree    = myfree   ;
      this->myzero    = myzero   ;
      this->growSize  = growSize ;
      this->blockSize = blockSize;
      this->pool_name = pool_name;
      this->error_message_out_of_memory = error_message_out_of_memory;
      this->error_message_cannot_grow   = error_message_cannot_grow  ;
      this->high_water_mark = 0;
      this->bytes_currently_allocated = 0;

      // Create the initial pool if the pool allocator is to be used
      pools.push_back( LinearAllocator( initialSize , blockSize , mymalloc , myfree , myzero ,
                                        pool_name , error_message_out_of_memory) );
    }


    /** @brief Finalize the pool allocator, deallocate all individual pools.
      */
    void finalize() {
      if (pools.size() > 0) {
        // SYCL rev6: The memory is freed without waiting for operating on
        // it to be completed. CUDA/HIP performs implicit synchronizations
        #ifdef YAKL_ARCH_SYCL
          fence();
        #endif
        if (yakl_mainproc()) std::cout << "Pool Memory High Water Mark:       " << get_high_water_mark() << std::endl;
        if (yakl_mainproc()) std::cout << "Pool Memory High Water Efficiency: " << get_pool_high_water_space_efficiency() << std::endl;
        if (! waiting_events.empty()) free_completed_waiting_entries();
        pools = std::list<LinearAllocator>();
        high_water_mark = 0;
        bytes_currently_allocated = 0;
      }
    }


    /** @brief **[USEFUL FOR DEBUGGING]** Print all allocations left in this pool object.
      * @details If you call this regularly in a loop, you can determine if an object fails to be deallocated
      *          properly by looking for repeated entries that grow in number each iteration. */
    void printAllocsLeft() {
      // Used for debugging mainly. Prints all existing allocations
      for (auto it = pools.begin() ; it != pools.end() ; it++) {
        it->printAllocsLeft();
      }
    }


    /** @brief Allocate the requested number of bytes using the requested label, and return the pointer to allocated space.
      * @details The pool allocator will search from beginning to end for a slot large enough to fit
      * this allocation request. This minimizes segmentation at the cost of a linear search time.
      * If the current pool(s) do not contain enough room, a new pool is created.
      *
      * Attempting to allocate zero bytes will return `nullptr`. This is a thread safe call.
      *
      * This always checks to see if entries waiting on stream events are able to be deallocated before allocating. */
    void * allocate(size_t bytes, char const * label="") {
      if (bytes == 0) return nullptr;
      if (! waiting_events.empty()) free_completed_waiting_entries();
      // Loop through the pools and see if there's room. If so, allocate in one of them
      bool room_found = false;  // Whether room exists for the allocation
      bool linear_bug = false;  // Whether there's an apparent bug in the LinearAllocator allocate() function
      void *ptr;                // Allocated pointer
      // Protect against multiple threads trying to allocate at the same time
      mtx1.lock();
      {
        // Start at the first pool, see if it has room.
        // If so, allocate in that pool and break. If not, try the next pool.
        for (auto it = pools.begin() ; it != pools.end() ; it++) {
          if (it->iGotRoom(bytes)) {
            ptr = it->allocate(bytes,label);
            room_found = true;
            if (ptr == nullptr) linear_bug = true;
            break;
          }
        }
        // If you've gone through all of the existing pools, and room hasn't been found, then it's time to add a new pool
        if (!room_found) {
          if (bytes > growSize) {
            std::cerr << "ERROR: For the pool allocator labeled \"" << pool_name << "\":" << std::endl;
            std::cerr << "ERROR: Trying to allocate " << bytes << " bytes (" << bytes/1024./1024./1024. << " GB), "
                      << "but the current pool is too small, and growSize is only "
                      << growSize << " bytes (" << growSize/1024./1024./1024. << " GB). \nThus, the allocation will never fit in pool memory.\n";
            std::cerr << "This can happen for a number of reasons. \nCheck the size of the variable being allocated in the "
                      << "line above and see if it's what you expected. \nIf it's absurdly large, then you might have tried "
                      << "to pass in a negative value for the size, or the size got corrupted somehow. \nNOTE: If you compiled "
                      << "for the wrong GPU artchitecture, it sometimes shows up here as well. \nIf the size of the variable "
                      << "is realistic, then you should increase the initial pool size and probably the grow size as "
                      << "well. \nWhen individual variables consume sizable percentages of a pool, memory gets fragmented, and "
                      << "the pool space isn't used efficiently. \nLarger pools will improve that. "
                      << "\nIn the extreme, you could create "
                      << "an initial pool that consumes most of the avialable device memory. \nIf that still doesn't work, then "
                      << "it sounds like you're choosing a problem size that's too large for the number of compute "
                      << "nodes you're using.\n";
            std::cerr << error_message_cannot_grow << std::endl;
            printAllocsLeft();
            die();
          }
          pools.push_back( LinearAllocator( growSize , blockSize , mymalloc , myfree , myzero ,
                                            pool_name , error_message_out_of_memory) );
          ptr = pools.back().allocate(bytes,label);
        }
        bytes_currently_allocated += ( (bytes-1)/blockSize+1 )*blockSize;
        high_water_mark = std::max( high_water_mark , bytes_currently_allocated );
      }
      mtx1.unlock();
      if (linear_bug) {
        std::cerr << "ERROR: For the pool allocator labeled \"" << pool_name << "\":" << std::endl;
        die("ERROR: It looks like you've found a bug in LinearAllocator. Please report this at github.com/mrnorman/YAKL");
      }
      if (ptr != nullptr) {
        return ptr;
      } else {
        std::cerr << "ERROR: For the pool allocator labeled \"" << pool_name << "\":" << std::endl;
        std::cerr << "Unable to allocate pointer. It looks like you might have run out of memory.";
        die( error_message_out_of_memory );
      }
      return nullptr;
    };


    bool thisIsMyPointer(void *ptr) {
      std::lock_guard lock(mtx1);
      // Go through each pool.
      for (auto it = pools.rbegin() ; it != pools.rend() ; it++) {
        if (it->thisIsMyPointer(ptr)) {
          return true;
        }
      }
      return false;
    }

    /** @brief Free the passed pointer, and return the pointer to allocated space.
      * @details Attempting to free a pointer not found in the list of pools will result in a thrown exception */
    void free(void *ptr , char const * label = "" ) {
      bool pointer_valid = false;
      // Protect against multiple threads trying to free at the same time
      mtx1.lock();
      {
        // Go through each pool. If the pointer lives in that pool, then free it.
        for (auto it = pools.rbegin() ; it != pools.rend() ; it++) {
          if (it->thisIsMyPointer(ptr)) {
            size_t bytes = it->free(ptr,label);
            bytes_currently_allocated -= bytes;
            pointer_valid = true;
            break;
          }
        }
      }
      mtx1.unlock();
      if (!pointer_valid) {
        std::cerr << "ERROR: For the pool allocator labeled \"" << pool_name << "\":" << std::endl;
        std::cerr << "ERROR: Trying to free an invalid pointer\n";
        die("This means you have either already freed the pointer, or its address has been corrupted somehow.");
      }
    };


    /** @brief Free the passed pointer, and return the pointer to allocated space.
      * @details Attempting to free a pointer not found in the list of pools will result in a thrown exception */
    void free_with_event_dependencies(void *ptr , std::vector<Event> events_in , char const * label = "") {
      mtx2.lock();
      waiting_entries.push_back( { std::string(label) , ptr , events_in } );
      for (int ind_event_in=0; ind_event_in < events_in.size(); ind_event_in++) {
        bool add_new_event = true;
        for (int ind_event_list=0; ind_event_list < waiting_events.size(); ind_event_list++) {
          if (events_in[ind_event_in] == waiting_events[ind_event_list]) add_new_event = false;
        }
        if (add_new_event) waiting_events.push_back( events_in[ind_event_in] );
      }
      mtx2.unlock();
    };


    /** @brief Check all deallcation entries that are waiting on stream events to see if those events have completed.
      *        If the events are completed, then free the entry from the pool. */
    void free_completed_waiting_entries() {
      mtx2.lock();
      // Loop through waiting events. Check if it's completed
      for (int ievent=0; ievent < waiting_events.size(); ievent++) {
        auto waiting_event = waiting_events[ievent];
        // If this event has completed, erase all instances of this event from all waiting entry event lists
        if (waiting_event.completed()) {
          auto &completed_event = waiting_event;
          // Loop through waiting entries and erase completed event from each entry's waiting event list
          for (int ientry=0; ientry < waiting_entries.size(); ientry++) {
            auto waiting_entry = waiting_entries[ientry];
            // Search through this waiting entry's waiting event list, and erase any occurrences of this completed event
            for (int i=0; i < waiting_entry.events.size(); i++) {
              if (waiting_entry.events[i] == completed_event) waiting_entry.events.erase(waiting_entry.events.begin()+i);
              break;
            }
            // if this waiting entry's waiting event list is empty, free its pointer & erase it from waiting entry list
            if (waiting_entry.events.empty()) {
              this->free( waiting_entry.ptr , waiting_entry.label.c_str() );
              waiting_entries.erase( waiting_entries.begin()+ientry );
            }
          }
          // Finally, erase this completed event from the waiting events list
          waiting_events.erase(waiting_events.begin()+ievent);
        }
      }
      mtx2.unlock();
    }


    /** @brief  Get the total capacity of all of the pools put together */
    size_t get_pool_capacity( ) const {
      size_t sz = 0;
      for (auto it = pools.begin() ; it != pools.end() ; it++) { sz += it->poolSize(); }
      return sz;
    }


    /** @brief Get the total number of allocations in all of the pools put together */
    size_t get_num_allocs( ) const {
      size_t allocs = 0;
      for (auto it = pools.begin() ; it != pools.end() ; it++) { allocs += it->numAllocs(); }
      return allocs;
    }


    /** @brief Get the current memory high water mark in bytes for all allocations passing through the pool */
    size_t get_high_water_mark() const { return high_water_mark; }


    /** @brief Get the current number of pools that have been allocated */
    int get_num_pools() const { return pools.size(); }


    /** @brief Get the current number of bytes that have been allocated in the pools (this is actual allocation, not pool capacity) */
    size_t get_bytes_currently_allocated() const {
      return bytes_currently_allocated;
    }


    /** @brief Get the current proportion of total capacity among pools that is actually allocated */
    double get_pool_space_efficiency() const {
      return static_cast<double>(get_bytes_currently_allocated()) / static_cast<double>(get_pool_capacity());
    }


    /** @brief Get the proportion of total capacity among pools that has been allocated at this largest past memory usage */
    double get_pool_high_water_space_efficiency() const {
      return static_cast<double>(get_high_water_mark()) / static_cast<double>(get_pool_capacity());
    }


  };

}
__YAKL_NAMESPACE_WRAPPER_END__

