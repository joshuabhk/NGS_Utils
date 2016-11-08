#!/usr/bin/env python

import sys
from sys import maxsize
from random import randint
from os.path import exists

class FASTQ_lane_unifier :
	'''
	This class unify the lanes to 1
	so that the illumina w
	'''
	def __init__( self, fn, outfn, simple=False ) :
		self.fn = fn
		self.tilexys = set()
		if exists(outfn) :
			print( "File %s already exits!"%outfn, file=sys.stderr )
			sys.exit(-1)
		else :
			outfp = open(outfn, 'w')

		self.count = 0
		if simple :
			self.run_simple( outfp )
		else :
			self.run( outfp )

	def run( self, outfp ) :
		for i, line in enumerate( open(self.fn) ) :
			if i%4 == 0 :
				outfp.write(self.uniquify_header(line))
			else :
				outfp.write(line)

	def run_simple( self, outfp ) :
		'''
		This is not the main method but wrote as a comparison purpose
		'''
		for i, line in enumerate( open(self.fn) ) :
			if i%4 == 0 :
				outfp.write(self.uniquify_header_simple(line))
			else :
				outfp.write(line)
				
	def uniquify_header( self, fastq_header, lane_number='1' ) :
		header1, header2 = fastq_header.split()
		l = header1.split(":")
		lane,tile,x,y = l[3], l[4], l[5],l[6]
		newid = ' '.join(l[4:7])
		if newid not in self.tilexys :
			self.tilexys.add(newid)
			l[3] = lane_number
			header1 = ':'.join(l)
			return '%s %s\n'%(header1, header2)
			
		oldid = newid
		while newid in self.tilexys :
			x, y = randint(1,maxsize), randint(1,maxsize)
			newid = '%s %s %s'%(tile,x,y)
		else :
			self.tilexys.add(newid)

		l[3] = lane_number	
		l[5] = str(x)
		l[6] = str(y)
		header1 = ':'.join(l)
		self.count += 1
		print( self.count, oldid, "==>>", newid )

		return '%s %s\n'%(header1, header2)


	def uniquify_header_simple( self, fastq_header, lane_number='1' ) :
		header1, header2 = fastq_header.split()
		l = header1.split(":")
		lane,tile,x,y = l[3], l[4], l[5],l[6]
		l[3] = lane_number
		header1 = ':'.join(l)
		return '%s %s\n'%(header1, header2)

if __name__=='__main__' :
	try :
		fq = sys.argv[1]
		ofq = sys.argv[2]

		FASTQ_lane_unifier( fq, ofq )
	except :
		print( "%s <Input> <Output>"%sys.argv[0] )
