#!/bin/bash
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash#answer-14203146

set -e

usage () {
    echo 'Combine Xena matrices by shared column (row name), i.e. grow horizontally.'
    echo 'usage: join_xena.sh [-h] [-o OUTPUT] file [file ...]'
    echo $'\npositional arguments:'
    echo ' file                   matrix(es) to be joined'
    echo $'\noptional arguments:'
    echo ' -h, --help             show this help message and exit'
    echo ' -o, --output OUTPUT    path to output file, including filename. Directory must'
    echo '                        exist and file must not exist (no overwritting).'
    exit 0
}

# Get the list of matrices to be joined
IFS=$'\n' # allows spaces in path
matrixlist=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -o|--output)
        if [ -e "$2" ]; then
            echo "Output file $2 exist! Overwrite is not supported."$'\n'
            usage
        fi
        outdir=$(dirname "$2")
        if [ ! -d "$outdir" ]; then
            echo "Output directory $outdir doesn't exist!"$'\n'
            usage
        fi
        output="$2"
        shift # past argument
        shift # past value
        ;;
        -h|--help)
            usage
        shift # past argument
        ;;
        *)    # unknown option; should be positional argument
        for path in "$1"
        do
            matrixlist+=("$path")
        done
        shift
        ;;
    esac
done

# Setup files
touch "$output"
temp="$output.$RANDOM.tmp"
sorted=$(mktemp -p "$outdir")
trap "{ rm -f $sorted; }" EXIT

# Real outer join one by one
for file in ${matrixlist[@]}
do
    echo "Merging $file ..."
    #LC_ALL=C sort -k 1b,1 "$file" > "$sorted"
    (head -n 1 "$file" && tail -n +2 "$file" | LC_ALL=C sort -k 1b,1) > "$sorted"
    LC_ALL=C join -t $'\t' -a1 -a2 -o auto --header --check-order "$output" "$sorted" > "$temp"
    mv "$temp" "$output"
done
