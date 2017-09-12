#!/usr/bin/env bash

wget http://ygob.ucd.ie/ygob/data/v7-Aug2012/Pillars.tab -O Pillars.tab
wget http://ygob.ucd.ie/ygob/data/v7-Aug2012/AA.fsa -O AA.fsa

mkdir -p queries
mkdir -p genomes

./ygob.py AA.fsa Pillars.tab
