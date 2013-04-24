BWT README
============
Author: Kemal Eren (kemal@kemaleren.com)

ABOUT
-----

`bwt` is a Python module for efficiently searching strings using the
Burrows-Wheeler transform.

Sample usage:

    import bwt
    import urllib2

    # get Pride and Prejudice, Chapters 1 through 5
    url = http://www.gutenberg.org/cache/epub/1342/pg1342.txt
    text = urllib2.urlopen(url).read()[677:30891]

    # pre-compute data structures
    bwt_data = bwt.make_all(text)

    # find all occurances of the word 'married', with up to three
    # mismatches.
    bwt.find('married', text, bwt_data=bwt_data)


Note: `bwt.make_all()` is not fast, because it uses a naive suffix
array algorithm. You can compute the suffix array seperately with any
efficient third-party module, and provide it to `make_all()`.


REQUIREMENTS
------------

* Python (tested with 2.7.4 and 3.3.1)
