Welcome to genephys' documentation
==================================

**genephys** is a Python tool to simulate electrophysiological data during a psychometric 
experiment, where for instance two stimuli are presented to a subject while
EEG or MEG data is recorded. The toolbox generates data both during baseline 
and after stimulation, continuously or in trials.
The data can be analysed using decoding analysis, and the results 
compared to decoding results from real data . 
To help with that, the package also provides basic decoding tools.

The generative model is described in `this paper <http://biorxiv.com>`_,
and a Jupyter notebook reproducing the main results is also available 
`here <https://github.com/vidaurre/genephys/tree/main/examples/jupyter_notebooks>`_.

Check :doc:`usage` section for further information.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   install
   rest
   task
   decoding
   graphics
   api
