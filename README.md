# Interference of Coherent Beams in Free Space: Separating Field Structure, Detector Response, and Back-Action

Reproducible code for the analytical results and figures of the article
*Interference of Coherent Beams in Free Space: Separating Field Structure, Detector
Response, and Back-Action* (Y. Y. Ghadi, S. Guizani, H. Hamam).

The article treats two mutually coherent beams crossing in free space and separates
three ingredients of any recorded fringe pattern: the prepared field, the observable
selected by the detector, and the back-action of the probe on the field. All results
are analytical; this repository verifies them numerically and regenerates every figure.
No experimental data are involved.

## Contents

| File | Purpose |
|---|---|
| `code/verify_benchmark.py` | Crossing-plane-wave benchmark: fields, Maxwell transversality, energy density, mean Poynting flux, preparation-mismatch bound (500 randomized checks; standard library only). |
| `code/crossed_design_check.py` | Crossed preparation × detector design: design-matrix rank, gauge deficit, χ² null distribution, detection of a non-factorizable interaction term, first-order expansion. |
| `code/identifiability_check.py` | Visibility model V = m κ (1 − 2w sin²θ) sinc(β sinθ) checked against direct averaging of partially coherent vector fields; coherence–contrast gauge; Jacobian rank; probe-size confounding; angle-dependent coherence; two-channel ratio test; least-squares covariance. |
| `code/verify_quantum_channel.py` | Selected two-mode absorption channel: single excitation shared by two modes, one excitation per mode, dark states. |
| `code/make_framework_overview.py` | Framework overview figure (PDF and SVG). |
| `code/make_figures.py` | Field profiles and signed fringe visibility versus crossing angle. |
| `code/make_amplitude_figure.py` | Addition of single-excitation detection amplitudes. |
| `code/make_identifiability_figure.py` | Angle-scan and two-channel ratio figure. |
| `Figures/` | Figures as produced by the scripts above. |

## Requirements and use

Python 3.9 or later with NumPy and Matplotlib (tested with Python 3.13, NumPy 2.5,
Matplotlib 3.11).

```
pip install -r requirements.txt
python code/verify_benchmark.py
python code/crossed_design_check.py
python code/identifiability_check.py
python code/verify_quantum_channel.py
python code/make_figures.py            # writes to Figures/
```

Each check script prints its checks and ends with a pass message; a failed check raises
an assertion error. All quantities are dimensionless (wavelength or amplitude set to 1).

## Citation

Please cite the archived release (see `CITATION.cff`; the DOI is assigned by Zenodo).

## License

MIT License (see `LICENSE`).
