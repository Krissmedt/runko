#!/bin/bash
#!/usr/bin/env python

echo "Running sdc_pypic.py w/ gyration_10.ini"
python3 sdc_pypic.py --conf gyration_10.ini

echo "Running sdc_pypic.py w/ gyration_20.ini"
python3 sdc_pypic.py --conf gyration_20.ini

echo "Running sdc_pypic.py w/ gyration_40.ini"
python3 sdc_pypic.py --conf gyration_40.ini

echo "Running sdc_pypic.py w/ gyration_80.ini"
python3 sdc_pypic.py --conf gyration_80.ini

echo "Running sdc_pypic.py w/ gyration_160.ini"
python3 sdc_pypic.py --conf gyration_160.ini

echo "Running sdc_pypic.py w/ gyration_320.ini"
python3 sdc_pypic.py --conf gyration_320.ini

echo "Running sdc_pypic.py w/ gyration_1000.ini"
python3 sdc_pypic.py --conf gyration_1000.ini
