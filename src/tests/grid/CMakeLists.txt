# ==========================================================================
list( APPEND grib_api_samples_data
reduced_gg_pl_160_grib1.tmpl  
reduced_gg_pl_512_grib2.tmpl    
regular_gg_sfc_grib2.tmpl
gg_sfc_grib1.tmpl                   
reduced_gg_pl_160_grib2.tmpl   
reduced_gg_pl_640_grib1.tmpl    
regular_ll_pl_grib1.tmpl
gg_sfc_grib2.tmpl                   
# reduced_gg_pl_2000_grib1.tmpl  
reduced_gg_pl_640_grib2.tmpl    
regular_ll_pl_grib2.tmpl
GRIB1.tmpl                          
# reduced_gg_pl_2000_grib2.tmpl  
reduced_gg_pl_80_grib1.tmpl     
regular_ll_sfc_grib1.tmpl
GRIB2.tmpl                          
reduced_gg_pl_200_grib1.tmpl   
reduced_gg_pl_80_grib2.tmpl     
regular_ll_sfc_grib2.tmpl
polar_stereographic_pl_grib1.tmpl   
reduced_gg_pl_200_grib2.tmpl   
reduced_gg_pl_grib1.tmpl        
rotated_ll_pl_grib1.tmpl
polar_stereographic_pl_grib2.tmpl   
reduced_gg_pl_256_grib1.tmpl   
reduced_gg_pl_grib2.tmpl        
rotated_ll_pl_grib2.tmpl
polar_stereographic_sfc_grib1.tmpl  
reduced_gg_pl_256_grib2.tmpl   
reduced_gg_sfc_grib1.tmpl       
rotated_ll_sfc_grib1.tmpl
polar_stereographic_sfc_grib2.tmpl  
reduced_gg_pl_320_grib1.tmpl   
reduced_gg_sfc_grib2.tmpl       
rotated_ll_sfc_grib2.tmpl
reduced_gg_ml_grib1.tmpl            
reduced_gg_pl_320_grib2.tmpl   
reduced_gg_sfc_jpeg_grib2.tmpl  
sh_ml_grib1.tmpl
reduced_gg_ml_grib2.tmpl            
reduced_gg_pl_32_grib1.tmpl    
reduced_ll_sfc_grib1.tmpl       
sh_ml_grib2.tmpl
reduced_gg_pl_1024_grib1.tmpl       
reduced_gg_pl_32_grib2.tmpl    
reduced_ll_sfc_grib2.tmpl       
sh_pl_grib1.tmpl
reduced_gg_pl_1024_grib2.tmpl       
reduced_gg_pl_400_grib1.tmpl   
regular_gg_ml_grib1.tmpl        
sh_pl_grib2.tmpl
reduced_gg_pl_1280_grib1.tmpl       
reduced_gg_pl_400_grib2.tmpl   
regular_gg_ml_grib2.tmpl        
sh_sfc_grib1.tmpl
reduced_gg_pl_1280_grib2.tmpl       
reduced_gg_pl_48_grib1.tmpl    
regular_gg_pl_grib1.tmpl        
sh_sfc_grib2.tmpl
reduced_gg_pl_128_grib1.tmpl        
reduced_gg_pl_48_grib2.tmpl    
regular_gg_pl_grib2.tmpl
reduced_gg_pl_128_grib2.tmpl        
reduced_gg_pl_512_grib1.tmpl   
regular_gg_sfc_grib1.tmpl
)

foreach( f ${grib_api_samples_data} )
  set(tgt data_${f})
  ecbuild_get_test_data( TARGET "${tgt}" NAME "${f}" DIRNAME ${PROJECT_NAME}/grid/grib_api/samples NOCHECK )
  list(APPEND tgts "${tgt}")
endforeach()

add_custom_target( atlas_grib_api_samples DEPENDS ${tgts} )

ecbuild_add_test( TARGET      atlas_test_GridSpec
                  BOOST
                  ARGS        reduced_gg_pl_160_grib1.tmpl
                  SOURCES     test_GridSpec.cc
                  CONDITION   ECKIT_HAVE_GRIB
                  LIBS        atlas 
                )


foreach( file ${grib_api_samples_data} )
    
    # Add the test which depends on the above executable, 
    # passing in individual grib samples files as ARGUMENTs
    ecbuild_add_test( 
                      TARGET test_grid_spec_${file}
                      TYPE SCRIPT
					  DEPENDS ${tgts}
                      COMMAND atlas_test_GridSpec
                      ARGS ${file}
                     )
endforeach()


# =========================================================================
                

#ecbuild_add_test(	TARGET      test_latlon
#					BOOST
#					SOURCES     test_latlon.cc
#					CONDITION   HAVE_GRID
#					LIBS        atlas  )

#ecbuild_add_test(   TARGET      test_gaussian
#					BOOST
#					SOURCES     test_gaussian.cc
#					CONDITION	HAVE_GRID
#					LIBS        atlas )

ecbuild_add_test(   TARGET      atlas_test_field
					BOOST
					SOURCES     test_field.cc
					CONDITION	ECKIT_HAVE_GRIB
					LIBS        atlas )

list( APPEND test_grid_from_grib_data
CMC_hrdps_east_TMP_TGL_2_ps2.5km_2014081512_P000-00.grib2
ll005005.grib
ll0101.grib
ll11.grib
ll55.grib
n1024.grib
#n1280.grib
n160.grib
n256.grib
n320.grib
n40.grib
n400.grib
n512.grib
n60.grib
rgg_n640.grib )

foreach( f ${test_grid_from_grib_data} )
  set(tgt data_${f})
  ecbuild_get_test_data( TARGET "${tgt}" NAME "${f}" DIRNAME ${PROJECT_NAME}/grid NOCHECK )
  list(APPEND tgts "${tgt}")
endforeach()

add_custom_target( atlas_tests_grid_get_gribs DEPENDS ${tgts} )

ecbuild_add_test(   TARGET      atlas_test_grid_from_grib
					BOOST
					SOURCES     test_grid_from_grib.cc
                    DEPENDS     atlas_tests_grid_get_gribs
					CONDITION	ECKIT_HAVE_GRIB AND OFF
					LIBS        atlas )
