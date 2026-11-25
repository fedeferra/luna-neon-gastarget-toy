# luna-neon-gastarget-toy

Straggling, resonant yield and detection efficiency in the LUNA ²²Ne extended gas target. Personal tool, not an official product of the LUNA collaboration.

Interactive, single page tool for the LUNA ²²Ne extended gas target (INFN LNGS). It shows, for a proton or alpha beam:

- energy and angular straggling along the target axis z;
- the resonant yield (Voigt profile folded with the beam energy distribution) versus beam energy;
- reactions produced and detected along z, using a Monte Carlo detection efficiency for three detectors (segmented NaI, BGO, HPGe).

Everything runs in the browser: no server, no dependencies. Open `index.html`, or use the online version
(GitHub Pages): https://fedeferra.github.io/luna-neon-gastarget-toy/.

## Use

Open `index.html` locally, or enable GitHub Pages (Settings, Pages, Source: GitHub Actions). The workflow in `.github/workflows/pages.yml` publishes the page at every push to `main`.

Main controls: beam energy E_b and the three helper buttons (E_b such that E = E_R at the target centre; E_b maximizing the reaction rate; E_b maximizing the detection rate), resonance energy and widths, projectile, detector dead time (detected events are multiplied by 1 − dead time), beam current, ²²Ne enrichment (100 % = pure ²²Ne; it scales only the absolute yield and rates), chamber length, gas stop position, beam heating, detector and photon energy.

## Model

- Density profile of the target: measured 2 mbar profile (EPJA 54 (2018) 44, Fig. 2), mirrored about z = 0 and scaled to 1 mbar; beam heating removes a fraction bh% of gas per 25 W of beam.
- Energy loss: proton stopping power in Ne; alpha obtained from the proton by Ziegler effective charge, normalized to 32.5e-15 eV cm² at 395 keV.
- Energy straggling: Bohr with a Lindhard Scharff type low velocity reduction, variance propagation d(σ²)/dn = Ω² − 2ε′σ². Angular straggling: Molière.
- Yield: dY/dz = n(z) · (enrichment) · Iσ · V(E(z) − E_R), V Voigt profile (Weideman Faddeeva function, N = 32), path integral with adaptive substeps.

## Photon lines

For the resonances with published branching ratios the detected rate for a selected line is the reaction rate × (emission probability per decay) × (peak efficiency at that energy and position). Choose one line (one photon at a time), a free photon energy (100 % per reaction), or your own list (`E_keV : intensity_% ; ...`).

Data: primary transitions from the ENSDF ²²Ne(p,γ)²³Na dataset (NNDC) for the resonances at 259.7, 436, 479, 629.5, 638.5, 725, 850.9, 896.6, 900.6, 919.8, 928.9, 947.8 and 980.9 keV (lab); the LUNA branching ratios of Depalo et al., PRC 94, 055804 (2016), for 156.2 and 189.5 keV. Secondary transitions are computed by following the cascade level by level with the gamma branching of the ²³Na adopted levels (NNDC; internal conversion neglected; levels without decay data, such as the 7385 keV one, are not followed). Intensities are photons per resonance decay (the sum exceeds 100 %); lines below 0.1 % are dropped. Each photon is treated as independent (no cascade summing). For ²²Ne(α,γ)²⁶Mg the 830 keV resonance has a list too (primary branching 54 % to the 1809 keV level and 46 % to the 7062 keV level from Hunt et al., PRC 99, 045804 (2019); cascade from the ²⁶Mg adopted levels, NNDC); the branching of Wolke et al. (1989) could not be retrieved. Resonances without published branching (upper limits, 661 keV, other alpha resonances) use the free energy or user defined lines. The level data were transcribed from the NNDC tables and should be checked against the source before quantitative use.

## Detection efficiency (Monte Carlo)

`mc/mc.c` (C, OpenMP) simulates photons emitted isotropically from points on the beam axis (21 positions, 0 to 30 cm; 28 energies, 20 keV to 20 MeV; 1e6 histories each). Photon cross sections come from xraylib up to 780 keV and Klein Nishina above, Bethe Heitler pair production (Geant4 style); no Rayleigh scattering, no photonuclear reactions, no escape of electrons or K X rays (bremsstrahlung as a single photon with the radiative fraction). The stainless steel beam tube (Ø20/24 mm, assumed) is included.

Peak efficiency: the deposit is smeared with a Gaussian of σ = FWHM(d)/2.3548 and counted in the window E ± FWHM(E). FWHM: NaI 7 % at 662 keV, BGO 10 %, HPGe 2 keV at 1332 keV, scaling as √E. Total efficiency: any deposit above 0.5 keV. The efficiency of the summed signal (addback) is used, without dead layers.

Geometries (cm), assumed where not stated by the experiment:

- NaI: six 10×10×20 crystals in a brick arrangement around a 10×10 square bore (z ± 10).
- BGO: six trapezoidal prisms (height 7) forming a hexagonal shell, bore apothem 3, outer apothem 10, length 28.
- HPGe: cylinder R = 4, h = 4, axis along x, centre at x = 8.

### Regenerating the efficiency table

```bash
cd mc
pip install numpy xraylib
python3 tables.py                      # writes tables.bin
gcc -O2 -fopenmp -o mc mc.c -lm
./mc 1000000 > eff.json                # about an hour on a few cores
cd ..
python3 tools/embed_sim.py mc/eff.json index.html
```

`mc/eff_smear.json` is the result embedded in the page.

## Limits

Indicative tool, not a substitute for a full Geant4 simulation: the geometries above are simplified, the Lindhard Scharff reduction is empirical, and the alpha stopping power is derived from the proton one.

## Licence

GPL-3.0, see `LICENSE`.
