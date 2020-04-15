### Eigen

ecbuild_find_package( NAME Eigen3 VERSION 3.3 QUIET )
ecbuild_add_option( FEATURE EIGEN
                    DESCRIPTION "Use Eigen linear algebra library"
                    CONDITION TARGET Eigen3::Eigen )

