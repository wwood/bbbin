#!/bin/bash

set -o pipefail

echo -n $1"	"
zcat $1 |
    awk 'NR%4==2{c++; l+=length($0)}
          END{
                print c"\t"l;
              }'
