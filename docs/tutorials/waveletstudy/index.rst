Tutorial: wavelet study
=======================

This example uses the preprocessed, epoched 250-Hz Flanker EEG on the
`Open Science Framework <https://osf.io/hgc65/>`_. Each input file should
contain one participant and one numeric congruency label per trial in
``EEG.etc.behavior.Cong``. Verify which two label values mean congruent and
incongruent in the study codebook. After ``nf_install``, put the following
code in a function and call it with ``nf_run(@my_analysis)``. Start in the
directory containing the preprocessed participant ``.set`` files.

.. code-block:: matlab

  files = dir('*.set');
  assert(~isempty(files),'No participant .set files found in this folder.');
  out = fullfile(pwd,'wavelet_results');
  if ~isfolder(out), mkdir(out); end
  useFiles = files(1:min(15,numel(files)));
  for i=1:numel(useFiles)
      EEG = pop_loadset('filename',useFiles(i).name,'filepath',useFiles(i).folder);
      EEG = nf_prepdata(EEG);
      TF = nf_tftransform(EEG,'method','dcwt','freqs',1:30,'times',[-0.4 1]);
      groups = [TF.behavior.Cong];
      assert(numel(groups)==TF.ntrls && numel(unique(groups))==2);
      TF = nf_normbase(TF,'TrialGroups',groups,'BaselineMethod','percent', ...
                       'BaselineTimes',[-0.4 -0.2]);
      TF = nf_resample(TF,-0.4:0.1:1,[]);
      [~,stem] = fileparts(useFiles(i).name);
      filename = fullfile(out,[stem '_tf.mat']);
      save(filename,'TF');
      processed(i) = dir(filename);
  end
  aggregate = nf_aggregate(processed,'aggtype','avg');
  nf_viewerapp(aggregate);

Check the label order in ``TF.normbase.conditions.values``. Averaged
aggregates store power as sensor × frequency × time × observation. Each
observation identifies its subject and condition. The following paired
comparison uses sensor 87, 5 Hz and 400 ms as in the original illustration.
Enter the condition codes verified in the codebook.

.. code-block:: matlab

  labels = TF.normbase.conditions.values;
  codes = input('Enter [congruent incongruent] numeric codes from the codebook: ');
  [present,indices] = ismember(codes,labels);
  assert(numel(indices)==2 && all(present) && indices(1)~=indices(2));
  [~,f] = min(abs(aggregate.freqs-5));
  [~,t] = min(abs(aggregate.times-0.4));
  paired = nan(aggregate.nsubjects,2);
  for subject=1:aggregate.nsubjects
      for c=1:2
          row = aggregate.trial_subjectindex(:)==subject & ...
                aggregate.trial_condindex(:)==indices(c);
          assert(nnz(row)==1);
          paired(subject,c) = aggregate.power(87,f,t,row);
      end
  end
  [h,p,ci,stats] = ttest(paired(:,2),paired(:,1));

The original 15-subject illustration reported t(14) = 2.97, p = .01,
95% CI [15.98, 99.22]. These historical plots identify a candidate effect;
choosing the sensor, time and frequency after viewing the same condition
contrast makes the follow-up test exploratory.

.. image:: nf_wavelet_flanker1.png
  :width: 400

.. image:: nf_wavelet_flanker3.png
  :width: 400

.. image:: nf_wavelet_flanker4.png
  :width: 400
