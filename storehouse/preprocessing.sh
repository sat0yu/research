#!/bin/sh

# pre-preprocessing
rm -rf \
    T1000S100.tsv T1000S1000.tsv \
    T10000S100.tsv T10000S1000.tsv T10000S10000.tsv \
    T100000S100.tsv T100000S1000.tsv T100000S10000.tsv

# create .tsv files
grep -e '^1000,100,' experiment_sliding.log \
    | tr ',' '\t' \
    > T1000S100.tsv

grep -e '^1000,1000,' experiment_sliding.log \
    | tr ',' '\t' \
    > T1000S1000.tsv

grep -e '^10000,100,' experiment_sliding.log \
    | tr ',' '\t' \
    > T10000S100.tsv

grep -e '^10000,1000,' experiment_sliding.log \
    | tr ',' '\t' \
    > T10000S1000.tsv

grep -e '^10000,10000,' experiment_sliding.log \
    | tr ',' '\t' \
    > T10000S10000.tsv

grep -e '^100000,100,' experiment_sliding.log \
    | tr ',' '\t' \
    > T100000S100.tsv

grep -e '^100000,1000,' experiment_sliding.log \
    | tr ',' '\t' \
    > T100000S1000.tsv

grep -e '^100000,10000,' experiment_sliding.log \
    | tr ',' '\t' \
    > T100000S10000.tsv

# create figures
gnuplot makeGraph.gnp
