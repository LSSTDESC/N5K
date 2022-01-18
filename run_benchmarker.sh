#!/bin/bash

mkdir -p outputs

python benchmarker.py confs/config_fkem.yml FKEM 200 > outputs/SNs_FKEM.txt
python benchmarker.py confs/config_matter.yml MATTER 200 > outputs/SNs_MATTER.txt
python benchmarker.py confs/config_levin.yml LEVIN 200 > outputs/SNs_Levin.txt
