#!/bin/bash

for file in ./logs/*.log; do 
	cat $file 
done
