#!/bin/bash
# wrapper script to enable code folding on travis-ci
#
# Note: Expects travis functions exported in before_script
#
# Usage: travis_wrap.sh fold-block-id script [script args]
# fold-block-id should not contain special characters, blanks newlines, ...

set -e
travis_time_finish
travis_fold start "$1"
    travis_time_start
        bash -c "${@:2:$#}"
    travis_time_finish
travis_fold end "$1"
