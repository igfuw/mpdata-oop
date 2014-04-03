#!/bin/bash

# generates test data for timing tests

Cx=.2
Cy=.5
it=3

for nxy in {0064,0081,0102,0128,0161,0203,0256,0323,0406,0512,0645,0813,1024,1290,1626,2048}; do
  nt=`echo "2^24/$nxy/$nxy"|bc`
  echo nx=ny=$nxy nt=$nt
  time ./cpp/inout-cpp $nxy $nxy $Cx $Cy $nt $it \
    1>"./tests/timing/nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:in" \
    2>"./tests/timing/nx=${nxy}_ny=${nxy}_Cx=${Cx}_Cy=${Cy}_nt=${nt}_it=${it}:out"
  if [ $? -ne 0 ]; then 
    echo "calling ./cpp/inout-cpp failed!"
    exit 1
  fi
done
