======
SPEPTIDE
======

speptide version 0.92

About
=====

Speptide is developed by a bioinformatics group from 
the SRI of Physical-Chemical Medicine.

Installation
======

Just run::

  $ make

in project directory

Usage
======

Run::

  $ speptide <exp mgf> <db mgf> <config>

Inputs::

  <exp mgf>    (string)  experiment spectra mgf file.
  <db mgf>     (string)  database spectra mgf file (with SEQ tag).
  <config>     (string)  ini file with params for search.

Output::

  [exp spectrum id] [db spectrum id] [position] [exp ami] [db ami] [db seq] [cos(theta)]

Configuration file
======

Configuration (**.ini**) file consists of a lot of different parameters. 
For detailed descriptions, see **default.ini** file and read doc/speptide.pdf documentation.

Example
======

Run::

  $ speptide tests/h1.mgf tests/h2.mgf params/default.ini 

Authors
=======

Dmitry Ischenko

Dmitry Alexeev

Contact
=======

For comments and requests, send an email to:

  ischenko.dmitry@gmail.com

