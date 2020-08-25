#!/bin/bash
#!/usr/bin/env python

echo "Running vv_pic.py"
python3 vv_pic.py --conf gyration.ini

echo "Running vv_pic_p.py"
python3 vv_pypic.py --conf gyration.ini
