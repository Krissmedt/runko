#!/bin/bash
#!/usr/bin/env python

echo "Running lf_pypic.py w/ helix_0025.ini"
python3 lf_pypic.py --conf helix_0025.ini

echo "Running lf_pypic.py w/ helix_0050.ini"
python3 lf_pypic.py --conf helix_0050.ini

echo "Running lf_pypic.py w/ helix_0100.ini"
python3 lf_pypic.py --conf helix_0100.ini

echo "Running lf_pypic.py w/ helix_0200.ini"
python3 lf_pypic.py --conf helix_0200.ini

echo "Running lf_pypic.py w/ helix_0400.ini"
python3 lf_pypic.py --conf helix_0400.ini

echo "Running lf_pypic.py w/ helix_0800.ini"
python3 lf_pypic.py --conf helix_0800.ini

echo "Running lf_pypic.py w/ helix_1000.ini"
python3 lf_pypic.py --conf helix_1000.ini

