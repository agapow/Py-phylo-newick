from setuptools import setup, find_packages
import os

version = '0.1'

setup(name='Py-phylo-newick',
      version=version,
      description="Routines for parsing and writing phylogenies in Newick format",
      long_description=open("README.txt").read() + "\n" +
                       open(os.path.join("docs", "HISTORY.txt")).read(),
      # Get more strings from
      # http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        "Programming Language :: Python",
        ],
      keywords='',
      author='',
      author_email='',
      url='http://svn.plone.org/svn/collective/',
      license='GPL',
      packages=find_packages(exclude=['ez_setup']),
      namespace_packages=['phylo'],
      include_package_data=True,
      zip_safe=False,
		test_suite='nose.collector',
      install_requires=[
          'setuptools',
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
)

