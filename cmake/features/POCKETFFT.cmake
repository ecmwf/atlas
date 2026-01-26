### POCKETFFT ...

if( atlas_HAVE_ATLAS_TRANS )

    set (_user_disabled_pocketfft OFF)
    if (DEFINED ATLAS_ENABLE_POCKETFFT AND NOT ATLAS_ENABLE_POCKETFFT)
        set(ENABLE_POCKETFFT ${ATLAS_ENABLE_POCKETFFT})
    endif()
    if (DEFINED ENABLE_POCKETFFT AND NOT ENABLE_POCKETFFT)
        set (_user_disabled_pocketfft ON)
    endif()

    if (NOT _user_disabled_pocketfft)
        # Find pocketfft header-only library!
        # If pocketfft is not found, the source codes will be downloaded, unless ENABLE_DOWNLOAD=OFF
        if (NOT DEFINED ENABLE_DOWNLOAD)
            set(ENABLE_DOWNLOAD ON)
        endif()
        if (NOT ENABLE_DOWNLOAD)
            ecbuild_find_package(pocketfft)
        else()
            ecbuild_find_package(pocketfft QUIET)
            if (pocketfft_FOUND AND EXISTS ${POCKETFFT_INCLUDE_DIRS})
                ecbuild_info("atlas FOUND pocketfft")
                ecbuild_info("   POCKETFFT_INCLUDE_DIRS : [${POCKETFFT_INCLUDE_DIRS}]")
            else()
                set(pocketfft_ROOT ${PROJECT_BINARY_DIR}/pocketfft)
                if (NOT EXISTS ${pocketfft_ROOT} )
                    set(pocketfft_COMMIT 0fa0ef591e38c2758e3184c6c23e497b9f732ffa) # from branch cpp, dated 30 Nov 2024
                    set(pocketfft_URL https://github.com/mreineck/pocketfft/archive/${pocketfft_COMMIT}.zip)
                    ecbuild_info("Downloading pocketfft (${pocketfft_URL}) ...")
                    file(DOWNLOAD ${pocketfft_URL} ${PROJECT_BINARY_DIR}/download/pocketfft.zip STATUS DOWNLOAD_STATUS)
                    # Separate the returned status code, and error message.
                    list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
                    list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)
                    # Check if download was successful.
                    if(${STATUS_CODE} EQUAL 0)
                        file(ARCHIVE_EXTRACT INPUT ${PROJECT_BINARY_DIR}/download/pocketfft.zip DESTINATION ${PROJECT_BINARY_DIR}/download)
                        file(REMOVE ${PROJECT_BINARY_DIR}/download/pocketfft.zip)
                        file(RENAME ${PROJECT_BINARY_DIR}/download/pocketfft-${pocketfft_COMMIT} ${pocketfft_ROOT})
                        ecbuild_info("Downloading pocketfft (${pocketfft_URL}) ... done")
                    else()
                        ecbuild_info("Downloading pocketfft (${pocketfft_URL}) ... FAILED: ${ERROR_MESSAGE}")
                    endif()
                endif()
                ecbuild_find_package(pocketfft)
            endif()
        endif()
    endif()

    ecbuild_add_option( FEATURE POCKETFFT
                        DEFAULT ON
                        DESCRIPTION "Enable use of pocketfft for FFT's used in Trans"
                        CONDITION pocketfft_FOUND )

endif()

if( NOT atlas_HAVE_POCKETFFT )
    unset( POCKETFFT_INCLUDE_DIRS )
endif()
