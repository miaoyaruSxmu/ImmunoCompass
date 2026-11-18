# ImmunoCompass 0.3.0

* Updated the package to the latest 781-sample Top390 panel reported in Table
  S6 and refitted the final SVR on the Table-S2-aligned 807-sample development
  set.
* Kept binary encoding as 1 for GeneA > GeneB and 0 otherwise, including ties.
* Added pair-specific development-mode encoding for unavailable gene pairs.
* Expanded `model_metadata` with panel provenance and the fitted SVR settings.

# ImmunoCompass 0.2.0

* Updated the bundled IBI model to the 807-sample development set defined in
  Table S2, including 26 GSE122459 samples.
* Added model-training provenance through `model_metadata`.
# ImmunoCompass 0.1.0

Initial package release with the locked 390-pair IBI model, example data, and
immune-context-aware BAS/CAS/IAS decomposition.
