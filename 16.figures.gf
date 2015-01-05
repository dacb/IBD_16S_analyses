#!/bin/bash

awk -F'\t' -v prefix="$0" -f $0.awk 13.summarize_persistent_otus.gf.percent.rehsaped.xls 
gnuplot $0.gplt

R --no-save < 16.figures.gf.MDS.R
R --no-save < 16.figures.gf.heatmap.R
