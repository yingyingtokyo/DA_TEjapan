#!/bin/bash
for i in $(seq -w 1 29); do
  for j in $(seq -w 0 23); do
    echo "${i}${j}" >> yr2019.txt
  done
done

