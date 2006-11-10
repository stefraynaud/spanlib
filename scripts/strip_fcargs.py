import re,sys

r = re.compile("^-l",re.I)
l = []
for lib in sys.argv[1:]:
	l.append(r.sub("",lib))

print str(l)
	
