#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Writing trees in Newick format.

"""

# TODO: is the string buffer really necessary?


__docformat__ = 'restructuredtext'


### IMPORTS

import re

from phylo.core import io

__all__ = [
	'Writer',
]


### CONSTANTS & DEFINES

_quotedNameRegex = re.compile (r"\s*'[^']*'\s*")
_spacesInNameRegex = re.compile (r"[\s\,]+")



### IMPLEMENTATION ####

class Writer (io.BaseWriter):
	"""
	A writer for Newick formatted phylogenies.

	"""
	# TODO: needs a prettyprint function?
	# TODO: need a 'just print structure' option?
	# TODO: a precision option for distances

	def __init__ (self, dialect=None):
		dialect = dialect or WriterDialect()
		io.BaseWriter.__init__ (self, dialect)

	def _write (self, in_tree, dest):
		"""
		Write out the passed phylotree object in Newick format.

		"""
		## Preparation:
		self._src_tree = in_tree
		self._dest_strm = dest
		self._dist_format = '%.3f'
		self._support_format = '%.2f'
		## Main:
		root = in_tree.root
		if (not root):
			root = in_tree.get_centroid_nodes()[0]
		self._writeNode (root)

	def _writeNode (self, node, parent=None):
		"""
		Write out the Newick format structure for this node.

		By definition, this includes all the nodes below.

		"""
		## Main:
		if (self._src_tree.is_node_tip (node)):
			# a simple (terminal) node
			name = node.get ('title') or node.get ('name')
			# if the name is not quoted and contains spaces, quote it
			if (not _quotedNameRegex.search (name)):
				if (_spacesInNameRegex.search (name)):
					name = "'%s'" % name
			self._dest_strm.write (name)
		else:
			# complex (internal) node
			self._dest_strm.write ('(')
			children = [n for n in self._src_tree.iter_adjacent_nodes(node) if
				(n is not parent)]
			first_node = True
			for child in children:
				if (first_node):
					first_node = False
				else:
					self._dest_strm.write (', ')
				self._writeNode (child, node)
			self._dest_strm.write (')')
			# do support value
			supval = node.get ('support', None)
			if (supval is not None):
				self._dest_strm.write (self._support_format % supval)
		# do the distance
		if parent:
			dist = self._src_tree.get_distance (node, parent)
		else:
			dist = node.get ('distance', None)
		if (dist is not None):
			self._dest_strm.write (':' + self._dist_format % dist)




### END ########################################################################


