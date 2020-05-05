#!/bin/bash
#!/usr/bin/env python

python3 vv_pic.py --conf gyration.ini
python3 stag_pic.py --conf gyration.ini
python3 pic.py --conf gyration.ini
