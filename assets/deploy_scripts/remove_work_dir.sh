#!/usr/bin/env bash
rsync --update -vrL results results_rsync2
rm -r work/
rm nextflow.nohup_2*
rm .nextflow.log*