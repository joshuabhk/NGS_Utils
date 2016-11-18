#!/usr/bin/env python

import os, sys
from pysam import Samfile

def filter_reads( sam, mapq, osam ) :
	for aread in sam :
		if aread.is_unmapped :
			continue

		if aread.mapq >= mapq :
			osam.write(aread)
			

if __name__ == '__main__' :
	mapq = int(sys.argv[1])
	for fn in sys.argv[2:] :
		isam = Samfile(fn)
		if fn.endswith( 'bam' ) :
			ofn = fn.replace( 'bam', 'map%s.bam'%mapq )
		elif fn.endswith( 'sam' ) :
			ofn = fn.replace( 'sam', 'map%s.bam'%mapq )
		else :
			ofn = fn + ".mapq%s.bam"%mapq

		if os.path.exists( ofn ) :
			print( "Error:", ofn, "already exists!" )
			continue

		osam = Samfile( ofn, 'wb', template=isam )
		filter_reads( isam, mapq, osam )
		osam.close()

