#!/usr/bin/env python

import sys
from  merge_lanes_uniqfy_tilexy import FASTQ_lane_unifier

if __name__=='__main__' :
	try :
		fq = sys.argv[1]
		ofq = sys.argv[2]

		FASTQ_lane_unifier( fq, ofq, simple=True )
	except :
		print( "%s <Input> <Output>"%sys.argv[0] )
