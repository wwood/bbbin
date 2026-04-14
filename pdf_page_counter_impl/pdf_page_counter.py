#!/usr/bin/env python3
import sys
from pypdf import PdfReader

for path in sys.argv[1:]:
    n = len(PdfReader(path).pages)
    print(f"{n}\t{path}")
