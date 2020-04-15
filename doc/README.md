# Building documentation

Required to be preinstalled:
- dot
- doxygen

##  m.css themed doxygen (Default, much preferred)

    module load python3  # 3.6 required
    cd <project_directory>
    virtualenv venv
    source venv/bin/activate
    pip install jinja2 Pygments
    git clone git://github.com/mosra/m.css
    export PATH=$(pwd)/m.css/documentation  # doxygen.py executable here

CMake option (not really required as this is default behaviour):

    -DATLAS_DOXYGEN_GENERATOR="m.css"

## Standard doxygen html

CMake option:

    -DATLAS_DOXYGEN_GENERATOR=stock

## Tweaking documentation:

#### Changing version displayed in documentation:

CMake option:

    -DATLAS_DOC_VERSION="latest"

#### Cache full path to doxygen or doxygen.py executable:

CMake option:

    -DATLAS_DOXYGEN_EXECUTABLE=<executable/path>

