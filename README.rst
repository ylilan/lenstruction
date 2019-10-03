
=============================
lenstruction- cluster lensing source reconstruction package
=============================
.. image:: https://badge.fury.io/py/cluster_lenstronomy.png
    :target: http://badge.fury.io/py/cluster_lenstronomy
.. image:: https://travis-ci.org/ylilan/cluster_lenstronomy.png?branch=master
    :target: https://travis-ci.org/ylilan/cluster_lenstronomy
``lenstruction``  is a package to perform source reonstruction in cluster strong lensing. For multiply-imaged source, the code correct lens model parameter  locally up to order flexion to improve the source reconstruction.  
The software is presented in Yang et al 2019 (in preparation). 
Installation
------------
git clone the repo https://github.com/ylilan/lenstruction.git 
.. code-block:: bash
    $cd /where/you/git/clone/lenstruction
.. code-block:: bash
    $python setup.py install --user
Requirements
------------
make sure you have already installed lenstronomy (>version 0.9.2), if not, please instal lenstronomy first.    
.. code-block:: bash
    $ pip install lenstronomy --user
Modelling Features
------------
Example notebooks
------------

Attribution
------------
If you make use of ``lenstruction``, please cite Yang et al 2019 (in preparation),
`Birrer & Amara 2018 <https://arxiv.org/abs/1803.09746v1>`_, and `Birrer et al 2015 <http://adsabs.harvard.edu/abs/2015ApJ...813..102B>`_

More  
------------
The code is open access, you can hack it or email me: <ylilan@astro.ucla.edu>.
