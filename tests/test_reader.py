#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for 'phylore.newick.reader', using nose.

"""


### IMPORTS ###

from StringIO import StringIO as Buffer 

from phylo.newick.reader import NewickReader


## CONSTANTS & DEFINES ###

BIG_TREE_FILE = "tests/in/big.tree"
BIG_TREE_2_FILE = "tests/in/big2.tree"
ANN_TREE_FILE = "tests/in/ann.tree"


### TESTS ###

class TestReadingBigTree (object):
	def setup(self):
		pass
	
	def test_bigread (self):
		rdr = NewickReader()
		t1 = rdr.read (BIG_TREE_FILE)
		t1._validate()
		t1._dump()
		assert (t1.count_nodes() == 182)
		assert (t1.count_branches() == 181)

	def test_big2 (self):
		rdr = NewickReader({'node_annotations': True})
		t1 = rdr.read (BIG_TREE_2_FILE)
		t1._validate()
		t1._dump()
	
	def teardown(self):
		pass



def test_annotated_tree():
	rdr = NewickReader({'node_annotations': True})
	t1 = rdr.read (ANN_TREE_FILE)
	t1._validate()
	t1._dump()
	assert (t1.count_nodes() == 182)
	assert (t1.count_branches() == 181)
	

### END ####################################################################