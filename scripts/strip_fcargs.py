import re,sys

r = re.compile("^-l",re.I)
l = []
for arg in sys.argv[1:]:
	for lib in re.split("\s+",arg):
		l.append(r.sub("",lib))

print str(l)
	
