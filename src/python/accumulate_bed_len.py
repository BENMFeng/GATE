#!/usr/bin/python
import sys
import types

k_dist = {}
if len(sys.argv) > 1:
	def usage():
		print >>sys.stderr, "python", sys.argv[0], "<IN:tab_file>"
	try:
		fp = open(sys.argv[1],"r")
	except IndexError:
		print >>sys.stderr, "No such file or directory"
		sys.exit(1)
	except IOError:
		usage()
		sys.exit(1)

else:
	fp=sys.stdin

for line in fp.readlines():
	col=line.split()
	if type(col[1])==type(col[2]):
		try:
			k_dist[col[0]] += int(col[2]) - int(col[1]) + 1
		except KeyError:
			k_dist[col[0]] = int(col[2]) - int(col[1]) + 1
fp.close()
for i in k_dist.keys():
	print i,k_dist[i]