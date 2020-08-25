#!/bin/bash
#!/usr/bin/env python

echo "Running vv_python_pic.py"
python3 vv_pypic.py --conf gyration.ini

echo "Running stag_python_pic.py"
python3 stag_pypic.py --conf gyration.ini
