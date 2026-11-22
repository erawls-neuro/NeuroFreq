nf_normbase
===========

Average power across trials and compute inter-trial phase coherence (ITPC)
from phase angles. For simple averaging with no baseline correction:

.. code-block:: matlab

  avg = nf_normbase(TF,'BaselineMethod','none');

To average by a numeric trial label and express each condition's power as
percent change from its own pre-event baseline:

.. code-block:: matlab

  groups = [TF.behavior.Cong];
  avg = nf_normbase(TF,'TrialGroups',groups,'BaselineMethod','percent', ...
                    'BaselineTimes',[-0.4 -0.2]);

Power is in ``avg.power``. When input angles exist, ``avg.itpc`` contains
the magnitude of the across-trial mean unit vector. The condition order
and contributing trial counts are in ``avg.normbase.conditions``. Supply
single-trial, nonnegative power for linear-scale inputs. Signed quadratic
distributions can be averaged directly along the trial dimension.
