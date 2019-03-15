#!/bin/bash
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash#answer-14203146

set -e

usage () {
    echo 'Combine Xena matrices by shared header row (1st row only; column name), i.e. grow vertically.'
    echo 'usage: join_xena.sh [-h] [-o OUTPUT] file [file ...]'
    echo $'\npositional arguments:'
    echo ' file                   matrix(es) to be joined'
    echo $'\noptional arguments:'
    echo ' -h, --help             show this help message and exit'
    echo ' -o, --output OUTPUT    path to output file, including filename. Directory must'
    echo '                        exist and file must not exist (no overwritting).'
    exit 0
}

colorder=`head -1 -q "$@" | awk '
BEGIN {
    ncol = 0
    FS = "\t"
    RS = "\n|\r\n"
}
{
    for (i = 1; i <= NF; i++) {
        if (! ($i in allfields)) {
            ncol++
            colorder[ncol] = $i
            allfields[$i] = ""
        }
    }
}
END {
    printf "%s", colorder[1]
    for (i = 2; i <= ncol; i++) {
        printf "\t%s", colorder[i]
    }
}
'`

awk -v c="$colorder" '
BEGIN {
    ncol = split(c, colorder, "\t")
    print c
    FS = "\t"
    RS = "\n|\r\n"
}
FNR==1 {
    for (i = 1; i <= NF; i++) {
        thisfields[i] = $i
    }
    next
}

{
    for (i in colorder) {
        output[colorder[i]] = ""
    }
    for (i = 1; i <= NF; i++) {
        output[thisfields[i]] = $i
    }
    printf "%s", output[colorder[1]]
    for (i = 2; i <= ncol; i++) {
        printf "\t%s", output[colorder[i]]
    }
    print ""
}
' "$@"
