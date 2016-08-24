#!/bin/bash
####################################################################################################
# RUN AUTOPEP8 ON ALL PYTHON FILES
# - runs autopep8 on all files with the .py extension
# - more information here: https://github.com/hhatto/autopep8
#                          https://www.python.org/dev/peps/pep-0008/
#
# USAGE:
# - run_autopep8.sh [optional source directory, defaults to the directory above scripts/]
####################################################################################################

#...................................................................................................
# CONFIGURATION
#...................................................................................................

DEFAULT_SOURCE_DIR=${PWD}/../
SOURCE_DIR=${DEFAULT_SOURCE_DIR:-1}
PYTHON_FILES=$(find ${SOURCE_DIR} -name "*py")
AUTOPEP8_OPTIONS=(--max-line-length 100
                  --verbose
                  --in-place)

#...................................................................................................
# RUN
#...................................................................................................

ABS_PATH=`cd ${SOURCE_DIR}; pwd`
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "RUNNING autopep8 ON ${ABS_PATH} WITH OPTIONS ${AUTOPEP8_OPTIONS[*]}"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
autopep8 ${AUTOPEP8_OPTIONS[*]} ${PYTHON_FILES}
