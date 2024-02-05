#!/bin/bash

set -e
#get the root dir (1st ancestor of the location where this script is stored)
SRC_DIR=$(dirname "$BASH_SOURCE")/..

# Invoke MIP on the rdlp instance files in paralel. All output is collected in a single output file
# Only output lines marked with a 'tag' are preserved in the final output. This is convenient to remove debug output.
function runBenchmark() {
    instances=($(fd --regex "wt(040|050|100).*[16]\.dat" .))
    machines=(4 2)
    settings=("settings/default_settings_8.json" "settings/default_settings_10.json" "settings/default_settings_12.json")
    printf "%s\n" "${machines[@]}" | xargs -I{} printf "%s {}\n" "${instances[@]}"| xargs -I{} printf "--json %s {}\n" "${settings[@]}" | parallel --no-notice -P 35 --eta --colsep ' ' "build/PM {}"
}

#switch to the root directory. This allows us to invoke this script from any directory. Then run the benchmark.
pushd $SRC_DIR
runBenchmark
popd
