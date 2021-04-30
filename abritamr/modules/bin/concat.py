#!/usr/bin/env python3
import pathlib, pandas, math, sys,  re


def aggregate_results(result_files, output):
    result_df = pandas.DataFrame()
    for r in result_files:
        tmp = pandas.read_csv(r, sep = '\t')
        if result_df.empty:
            result_df = tmp
        else:
            result_df = result_df.append(tmp)
    result_df.to_csv(output, sep = '\t')


result_files = []

for i in sys.argv[2:]:
    l = i.strip('[,]')
    # print(l)
    result_files.append(l)

aggregate_results(result_files = result_files, output = sys.argv[1])