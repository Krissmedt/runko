#!/bin/bash
#!/usr/bin/env python

echo "Running vv_pic.py"
python3 vv_pic.py --conf gyration.ini

echo "Running stag_pic.py"
python3 stag_pic.py --conf gyration.ini
