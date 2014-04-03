#!/bin/bash

cat Testing/Temporary/LastTest.log | awk -f filter0.awk \
  | egrep "^time:| of all | insns per cycle " \
  | awk -f filter1.awk | tr -s " " | tr -s "%" | tr "%" "0" \
  | awk -f filter2.awk \
  | awk -f filter3.awk 

