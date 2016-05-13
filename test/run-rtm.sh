#!/bin/bash
set -x

< hom.rsf ../rtm trace=tracediff.rsf isz=3 isxbeg=200 isxend=200 iskip=20 igz=2 igxbeg=101 igxend=301 --readwrite=y > rtm.rsf
