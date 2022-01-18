#!/bin/bash

mkdir -p outputs

python timer.py confs/config_fkem.yml FKEM > outputs/times_FKEM.txt
python timer.py confs/config_matter.yml MATTER > outputs/times_MATTER.txt
python timer.py confs/config_levin.yml LEVIN > outputs/times_Levin.txt
