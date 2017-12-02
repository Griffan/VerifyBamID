#!/bin/bash
set -e

# Default parameter
input=input.txt
output=output.txt
refType='1000g'
isGrey='TRUE'

# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
-------------------------------------------------------------------------------
Filename:    run.plot.sh
Revision:    1.0
Date:        2017/6/24
Author:      Fan Zhang
Email:       fanzhang@umich.edu
Website:     https://griffan.github.io/VerifyBamID/
Description: This script is for PC coordinates plot
-------------------------------------------------------------------------------
# 1. Input.txt, file contains PC coordinates of target sample
For example:
SampleID	PC1	PC2
WT.1	0.001	0.001
WT.2	-0.003	-0.003
WT.3	0.004	0.034

# 2. Output file, prefix of plot pdf file
TargetSample.pdf

OPTIONS:
    -r use 1000g or hgdp as reference panel samples, default:1000g
    -g plot reference panel sample as grey points, default:TRUE
    -i input file, file contains PC coordinates of target sample
    -o output file or output directory, default:output
    -h/? show help of script
Example:
    run.plot.sh -i input.txt -o TargetSample -r 1000g -g grey
EOF
}
if [ $# -eq 0 ]
then
	echo "No arguments supplied"
	usage
	exit 1
fi
# Analysis parameter
while getopts "r:g:h:i:o:" OPTION;do
  case $OPTION in
    r)
      refType=$OPTARG
      ;;
    h)
      usage
      exit 1
      ;;
    i)
      input=$OPTARG
      ;;
    o)
      output=$OPTARG
      ;;
    g)
      isGrey=$OPTARG
      ;;
    ?)
      usage
      exit 1
      ;;
  esac
done
my_dir="$(dirname "$0")"
Rscript $my_dir/Plot.R $input $output $refType $isGrey
