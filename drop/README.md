# Droplet cases

This directory groups the repository's droplet simulations by physical setup:

- [`pinned/`](pinned/README.md) — deterministic contact-angle/equilibrium cases with a precursor film and disjoining pressure. “Pinned” is the project label for this contact-angle study; it is not a newly imposed fixed-contact-line boundary condition.
- [`droplet-spreading/`](droplet-spreading/README.md) — free-spreading cos⁴ droplet solver. It is deterministic by default and includes opt-in white-noise and OU-SE modes.
- [`droplet-coalescence/`](droplet-coalescence/README.md) — fixed-AMR double-drop experiment with a history-dependent noise switch and a notebook for bridge height, centerline sections and top views.

The cases are kept separate because their questions, source implementations, and validation evidence differ.
