# MDSplus node reference for EAST

Consolidated mapping of diagnostic quantities to MDSplus nodes (global plasma
parameters, EFIT, ECE, XCS, radiation, divertor, gas injection and auxiliary
heating).

- **Full manual**: see [`mdsplus_extraction_manual.md`](mdsplus_extraction_manual.md)
  for the complete extraction workflow, node tables and Python scripts used to
  download the time traces from the ASIPP MDSplus server.
- **Formal reference**: Appendix B of the thesis.

## Access notes

The MDSplus server (`202.127.204.12`) is only reachable from the ASIPP local
network in Hefei. Access from outside is arranged through a remote-desktop
session to an ASIPP workstation. Contact the ASIPP impurity spectroscopy group
to coordinate access.

## Verified shot range

All node names in the manual were verified experimentally against shots
`#143064`–`#143079` (used in the Liu 2026 density-free regime study).
Node availability may vary in later campaigns.
