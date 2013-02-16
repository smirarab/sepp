import urllib2
import Bio.Phylo as bp
import os


# download the taxonomy archive
url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
filename = url.split('/')[-1]
if os.path.exists(filename):
    print 'Using existing copy of %s' % filename
else:
    r = urllib2.urlopen(urllib2.Request(url))
    assert r.geturl() == url
    with open(filename, 'wb') as output_file:
        output_file.write(r.read())
    r.close()

