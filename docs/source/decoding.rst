Decoding analysis
=================

After sampling, the goal is to put the data through the lens of decoding analysis
in order to compare with real data. In particular, we want to look at 
the temporal generalisation matrix (TGM)
a (``T`` by ``T``) matrix of decoding accuracies,
such that one decoding model is trained per time point and tested on each one of 
the time points of the trial in a cross-validation fashion. 
The diagonal of the TGM reflects how well can we can decode information time point by time point, 
which is often interpreted as 
reflecting the waxing and waning of the different stages of stimulus processing (or decision making). 
The off-diagonal of the TGM shows how well decoding models generalize to time points 
different to those where they were trained; and, therefore, 
is argued to reflect the stability of the neural code for the neural representation under study. 

Although the generative model described above can only sample discrete stimuli,
the decoders here can deal with both discrete and continuous data, i.e. classification and regression.
Of course, they can be applied to real data just as well, 
as far as they are in the same *time by trials by channels* format.
For regression, linear discriminant analysis is used;
for classification, ridge regression is used. 
Importantly, the decoders assume that only one stimulus is presented per trial. 

We first create a ``Decoder`` object, 
where we specify if the decoding problem is a classification one, as well as other options. 

.. code-block:: console

    decoder = decoders.Decoder(classification=True,binsize=1)

The options for the initialisation of the ``Decoder`` object are

* ``classification``: is this a classification problem? 
  Default: True. 
* ``get_TGM``: whether we get a (``T`` x ``T``) TGM or just a vector of accuracies with ``T`` elements
  (i.e. without generalising). 
  Default value: ``True``.
* ``binsize``: the size of the bin or window on which to perform the training and the predictions.
  Default value: 1 (that is, time point by time point).
* ``binstep``: by how many time points we slide the window, 
  if we are doing a sliding window prediction (by specifying ``binsize > 1``).
  Default value: 1.
* ``cvscheme``: the (optional) cross-validation scheme. 
  This is logical numpy array of dimensions (no. of trials by no. of cross-validation folds).
  If an element ``[i,j]`` is ``True``, that means that the i-th trial will be using for testing
  in the j-th fold; if it is ``False``, it will be used for training. 
* ``ncv``: number of cross-validation folds (only used if ``cvscheme`` was not supplied).
  Default value: 10.
* ``alpha``: regularisation parameter. Default: 0.01. 

For example, we could manually create a cross-validation scheme as 

.. code-block:: console

    test_inds = np.zeros((N,10),dtype=bool)
    for icv in range(10):
        testing = (np.arange(N/10) + icv * (N/10)).astype(int)
        test_inds[testing,icv] = True
    cvscheme = test_inds

And specify it as an option in the initialisation of the object.
This might be useful if there is any specific structure in the trial set we wish to account for. 

Then, we can run the decoder: 

.. code-block:: console

    decoder = decoders.Decoder(classification=True,cvscheme=cvscheme)
    accuracy,betas = decoder.decode(X,stim)

The output has two variables, ``accuracy`` and ``betas``. 

The first is obviously the accuracy of the cross-validated decoder.
If we are computing a TGM (option ``get_TGM=True``)
``accuracy`` will be 
(``T`` by ``T`` by ``(Q*(Q-1)/2)``)-dimensional,
if this is a classification problem with more than two conditions;
or (``T`` by ``T``)-dimensional 
if this is a two-condition classification problem or a regression problem (a continuous stimulus).
If ``get_TGM=False``, then ``accuracy`` is either (``T`` by ``(Q*(Q-1)/2)``)
or a 1D array with ``T`` elements, for the ``Q>2`` case and otherwise. 
For classification, accuracy is the proportion of correctly predicted classes;
for regression, this the correlation between the ground-truth and the predicted stimulus values.

The second output argument is ``betas``, containing the decoding coefficients per time point or bin.
Its dimension is *channels* by *bins/time points* by ``(Q*(Q-1)/2)``,
or just *channels* by *bins/time points* if this is two-class classification or regression.