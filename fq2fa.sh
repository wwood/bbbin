#!/bin/bash
awk '{print ">" substr($0,2);getline;print;getline;getline}'
