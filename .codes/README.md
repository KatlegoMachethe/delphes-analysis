# Delphes analysis

Build from this directory with `make`. The executable is written to
`/home/katlego/wits/phd/ttaxion/bin/analyzedelphes` by default.
The build uses `ROOTSYS/bin/root-config` (or `root-config` on `PATH`) and
the Delphes installation at `~/wits/phd/software/Delphes-3.5.1`.
For another installation, pass `ROOT_CONFIG=/path/to/root-config` and
`DELPHES_DIR=/path/to/Delphes` to `make`.

Run from the parent data directory:

```bash
cd ..
../bin/analyzedelphes codes/params.json
```

Alternatively, `make run` from `codes` builds the executable and runs it
from the parent directory. With no argument, the executable reads
`params.json` in its current working directory.

`params.json` lists the 10, 50 and 100 GeV signal files in that order.
Relative file names are resolved against `base_dir`. The current configuration
uses `sig_10GeV.root`, `sig_50GeV.root` and `sig_100GeV.root`, with background
enabled from `bg.root`. The `delphes_sig_*.root` files in this dataset have no events with
two reconstructed photons and therefore yield empty histograms for this analysis.

The analysis reads nominal weights from `Event` as `HepMCEvent` and requires
two photons, at least one jet with `BTag == 1`, and at least two light jets.
It reports skipped and selected event counts. Only the hadronic final state
is implemented.

The reconstructed top is the sum of the leading b-tagged jet and the two
leading light jets. Only the b-jet receives the exact sequence from
`/home/katlego/wits/phd/ttaxion/analysis/newcode/histo.cpp`:

```cpp
bbar_v.RotateZ(-tbar_v.Phi());
bbar_v.RotateY(-tbar_v.Theta());
bbar_v.Boost(-tbar_v.BoostVector());
```

The top is unrotated. Photons and light jets keep their lab four-vectors.
`cost_bbar` uses the transformed b-jet's `CosTheta()`, while the three
azimuthal differences involving it compare the lab photon/jet with that
transformed b-jet, exactly as in the reference. A star on the b-jet plot
labels denotes this reference transformation. The other three observables
and all photon selection cuts continue to use lab-frame quantities.

Seven signal text histograms, each containing a header and 20 bins, are
written to each of `../histograms/m10`, `../histograms/m50` and
`../histograms/m100`. The background is processed once and produces seven
shared files in `../histograms/background`, named `<observable>_background.txt`.
A rerun replaces those files.
To save the analysis log from `codes`, use `make run > ../histograms/run.log 2>&1`.

The existing selection cuts are retained, including the **57–63 GeV** mass
window for the 50 GeV sample, as requested. Azimuthal differences retain the
original subtraction without wrapping: values outside the histogram range
enter underflow/overflow bins, whose contents are printed but excluded from
the text tables. Bin contents are weighted yields, not normalized distributions.
If `background` is enabled, it uses the common selection only: two photons,
at least one b-tagged jet, at least two light jets, and both photon eta values
between 0.5 and 5.0. No mass-dependent photon-PT thresholds or diphoton mass
windows are applied to the background. Its weights use its own input event
count, and the same distribution is shared by all three signal masses.
Signal and background tables retain weighted bin contents; no unit-area
normalization is written to these files.

Run `python3 ../plot_histograms.py` from `codes` to regenerate the PNGs.
The plotting script normalizes each curve to unit area over the saved bins.
Each per-mass plot overlays the same shared background for that observable.
The combined plots in `../histograms/comparisons` display the three signals
and exactly one dotted black background curve. Background-only PNGs are
also saved alongside the shared background text files.

Run `make test` to check the reference sequence against an analytic off-axis
example, verify unchanged input vectors and invariant masses, and test
invalid-top rejection.
