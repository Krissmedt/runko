#!/bin/bash
#!/usr/bin/env python

echo "Running vv_pypic.py w/ penning_0050.ini"
python3 vv_pypic.py --conf penning_0050.ini

echo "Running vv_pypic.py w/ penning_0100.ini"
python3 vv_pypic.py --conf penning_0100.ini

echo "Running vv_pypic.py w/ penning_0200.ini"
python3 vv_pypic.py --conf penning_0200.ini

echo "Running vv_pypic.py w/ penning_0400.ini"
python3 vv_pypic.py --conf penning_0400.ini

echo "Running vv_pypic.py w/ penning_0800.ini"
python3 vv_pypic.py --conf penning_0800.ini

echo "Running vv_pypic.py w/ penning_1600.ini"
python3 vv_pypic.py --conf penning_1600.ini

