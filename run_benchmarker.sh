#!/bin/bash

mkdir -p outputs

echo "FKEM"
python benchmarker.py confs/config_fkem.yml FKEM 200 none
python benchmarker.py confs/config_fkem.yml FKEM 200 full
python benchmarker.py confs/config_fkem.yml FKEM 200 half
python benchmarker.py confs/config_fkem.yml FKEM 200 quarter

echo "MATTER"
python benchmarker.py confs/config_matter.yml MATTER 200 none
python benchmarker.py confs/config_matter.yml MATTER 200 full
python benchmarker.py confs/config_matter.yml MATTER 200 half
python benchmarker.py confs/config_matter.yml MATTER 200 quarter

echo "Levin"
python benchmarker.py confs/config_levin.yml LEVIN 200 none
python benchmarker.py confs/config_levin.yml LEVIN 200 full
python benchmarker.py confs/config_levin.yml LEVIN 200 half
python benchmarker.py confs/config_levin.yml LEVIN 200 quarter
