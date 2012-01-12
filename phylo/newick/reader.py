#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Class for reading trees in Newick format.

"""
# TODO: parse scientific format, e.g. 9.58443e-05

__docformat__ = 'restructuredtext'


### IMPORTS

import re, cStringIO

from phylo.core import io
from phylo.core.tree import Tree


### CONSTANTS & DEFINES ###

_quotedNameRegex = re.compile (r"\s*'([^']*)'\s*")
_nameRegex = re.compile (r"\s*([a-zA-Z0-9\-_\?\*\/]+)\s*")
_cleanSpaceRegex = re.compile (r"\s+")
_distRegex = re.compile (r"\s*:\s*(\-?\d\.\d+[eE]\-\d+|\-?0?\.\d+|\d+\.\d+|\d+)\s*",
	re.IGNORECASE)
_SUPPORTVAL_RE = re.compile (r"\s*(\-?0?\.\d+|\d+\.\d+|\d+)\s*")
_NODE_ANN_REGEX = re.compile (r"^\s*\[\&([^\]]*)\]\s*")
_ANN_LIST_SPLIT_RE = re.compile (r',([a-zA-Z]+)')



### IMPLEMENTATION ###

class NewickReader (io.BaseReader):
	"""
	A parser for Newick formatted phylogenies.

	Some interpretation of the Newick format is required. This parser allows:

	* singular / singleton nodes, e.g. `((a))`

	* un-named tips, e.g. `(,(,),)`

	* the naming of internal nodes, e.g. `((ab,cd)ef:0.5)gh:0.4`

	* Names consisting of alphanumeric characters and underscores or those that
	  are within single quotes.

	Things we don't precisely do according to the standard:

	* "Unquoted labels may not contain blanks, parentheses, square brackets,
	  single_quotes, colons, semicolons, or commas." (What exactly is left?)

	* "Underscore characters in unquoted labels are converted to blanks."
	  (Interpreted as blanks, I take it.)

	* "Single quote characters in a quoted label are represented by two single
	  quotes." (Too hard.)

	* "Blanks or tabs may appear anywhere except within unquoted labels or
	  branch_lengths." (We actually collapse all consecutive whitespace to a
	  single space.)

	* "Comments are enclosed in square brackets and may appear anywhere
	  newlines are permitted.". (No. Too hard.)

	For example::
	
		>>> from StringIO import StringIO
		>>> s1 = "(A, (B, C))"
		>>> rdr = NewickReader()
		>>> t1 = rdr.read (StringIO (s1))
		>>> t1.count_nodes(), t1.count_branches()
		(5, 4)
		>>> s2 = "('Foo bar':1.11, (B, C):2.76e-1):0.123"
		>>> t2 = NewickReader().read (StringIO (s2))
		>>> t2.root['distance']
		0.123
		
		
	"""
	# TODO: should be a more efficient way of parsing.
	# TODO: should pass in some parsing options (whitespace treatment etc.)
	# TODO: it would also be nice later to have this so it can read multiple
	# trees and behaves like a proper parse, i.e. like the CSV reader:
	#   theReader = csv.reader (theCanonStr)
	#   for line in theReader:
	#      etc.
	# TODO: It would be nice to get rid of all internal space, but what if
	# there are quoted names?

	def __init__ (self, dialect={}, tree_kls=Tree):
		io.BaseReader.__init__ (self, dialect=dialect)
		self._tree_kls = tree_kls

	def _read (self, src):
		"""
		Read the passed tree and return a Phylotree structure.

		:Params:
			src
				An open file or file-like object.

		:Returns:
			A phylotree.

		"""
		## Preconditions & preparation:
		tree_str = src.read()
		# clean up string
		# get rid of flanking whitespace and trailing semi-colon, if any
		tree_str = tree_str.strip()
		if (tree_str.endswith (';')):
			tree_str = tree_str[:-1].strip()
		# compact excessive internal space
		tree_str = _cleanSpaceRegex.sub (' ', tree_str)
		# validate tree string
		assert (tree_str), "null string passed as tree data"
		theNumLeftBraces = tree_str.count ("(")
		theNumRightBraces = tree_str.count ("(")
		assert (theNumLeftBraces == theNumRightBraces), \
			"the number of left and right braces (%s, %s) don't match" % (
			theNumLeftBraces, theNumRightBraces)
		# NOTE: we can't make any assertion about the number of commas, because
		# (A,B,C,D,E) and ((((A,B)))) are both legal trees.
		#theNumCommas = tree_str.count (",")

		## Main:
		self._curr_src = cStringIO.StringIO (tree_str)
		self._curr_tree = self._tree_kls()
		self._parse_node (None)

		## Postconditions & return:
		return self._curr_tree

	def _peek (self, look_ahead=1):
		# TODO: need a peek char and peek line?
		theLookAhead = self._curr_src.read (look_ahead)
		self._curr_src.seek (-look_ahead, 1)
		return theLookAhead

	def _consumeLeadingSpace (self):
		try:
			while (self._curr_src.read(1) == ' '):
				pass
			self._curr_src.seek (-1, 1)
		except StopIteration:
			return

	def _parseName (self):
		"""
		Returns the on the name on the front of the string, if it's there.
		"""
		theNextStr = self._curr_src.getvalue()[self._curr_src.tell():]
		theMatch = _quotedNameRegex.match (theNextStr) or \
			_nameRegex.match (theNextStr)
		if (theMatch):
			self._curr_src.read (len (theMatch.group (0)))
			return theMatch.group(1)
		else:
			return None

	def _parseDist (self):
		theNextStr = self._curr_src.getvalue()[self._curr_src.tell():]
		theMatch = _distRegex.match (theNextStr)
		if (theMatch):
			self._curr_src.read (len (theMatch.group (0)))
			return float (theMatch.group(1))
		else:
			return None
			
	def _parse_support_value (self):
		theNextStr = self._curr_src.getvalue()[self._curr_src.tell():]
		theMatch = _SUPPORTVAL_RE.match (theNextStr)
		if (theMatch):
			self._curr_src.read (len (theMatch.group (0)))
			return float (theMatch.group(1))
		else:
			return None

	def _parse_node (self, par_node):
		"""
		Parse the current node in the passed string into the growing tree.

		Note that we allow singleton nodes - nodes with a single child - in this
		implementation, ie. `((a))`.
		"""
		## Main:
		# create new node
		if (par_node is None):
			new_node = self._curr_tree.add_root()
		else:
			new_node, new_branch= self._curr_tree.add_node (par_node)
		
		# move to start of tree string
		self._consumeLeadingSpace()
		if (self._peek() == '('):
			# if it's a compound node
			# consume leading brace
			self._curr_src.read(1)
			# parse the first child node
			self._parse_node (new_node)
			# now we look for one or more following nodes
			self._consumeLeadingSpace()
			#print self._curr_tree._nodes
			while (self._peek() == ','):
				self._curr_src.read(1)
				self._parse_node (new_node)
				self._consumeLeadingSpace()
			# consume end bracket
			self._consumeLeadingSpace()
			posn = self._curr_src.tell()
			assert (self._peek() == ')'), \
				"can't find end brace in '%s...' at position %s '%s'" % (
					self._curr_src.getvalue ()[:10],
					posn,
					self._curr_src.getvalue ()[posn:posn+10],
				)
			self._curr_src.read (1)
		# read support val if it's there
		# TODO: this eats the name if it's numeric
		#support_val = self._parse_support_value()
		#if (support_val is not None):
		#	new_node['support'] = support_val
		# read node name
		node_name = self._parseName()
		if (node_name is not None):
			new_node.title = node_name
		#print self._curr_src.getvalue()[:30]
		if self.dialect.get ('node_annotations', True):
			x = self._parse_node_annotations (new_node)
		# read distance (if it's there)
		dist_to_node = self._parseDist()
		if (dist_to_node is not None): 
			if par_node:
				new_branch['distance'] = dist_to_node
			else:
				new_node['distance'] = dist_to_node
		return new_node
	
	def _parse_node_annotations (self, curr_node):
		theNextStr = self._curr_src.getvalue()[self._curr_src.tell():]
		theMatch = _NODE_ANN_REGEX.match (theNextStr)
		if (theMatch):
			self._curr_src.read (len (theMatch.group (0)))
			# TODO: must be a better way to split this string (has internal commas)
			anns = _ANN_LIST_SPLIT_RE.sub (r'||\1', theMatch.group(1)).split ('||')
			prs = dict ([x.split ('=', 1) for x in anns])
			curr_node.update (prs)
		else:
			return None		

### TEST & DEBUG ###

def _doctest ():
   import doctest
   doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
   _doctest()


### END ######################################################################
