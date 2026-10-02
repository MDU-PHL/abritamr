# abriTAMR V2

abriTAMR V2 provides separate commands for detecting AMR mechanisms, assigning
AMR status, and creating reports. The main command groups are:

- [`scan`](commands/scan.md): run AMRFinderPlus on assemblies or annotate existing
  AMRFinderPlus output.
- [`amr_status`](commands/amr_status.md): add abriTAMR reporting status and AMR
  typing criteria to scan results.
- [`linelist`](commands/linelist.md): create a sample-level AMR report.
- [`matrix`](commands/matrix.md): create a presence/absence matrix.
- [`infer`](commands/infer.md): infer genotypic drug susceptibility (gDST) where
  species-specific rules are available.
- [`run`](workflow.md): run these steps together as a complete workflow.
- [Catalog and rules](catalog-and-rules.md): generate a reference gene catalog
  and gDST rules.

Run `abritamr --help` or `abritamr <command> --help` for the options supported
by your installed version. Commands write CSV to standard output by default;
use `--format tab` for tab-delimited output, and shell redirection to save it.

For assembly inputs, AMRFinderPlus and the abriTAMR dependencies must be
available in the environment. Supply a species with `--species` when known;
otherwise abriTAMR attempts to identify it for scan and downstream inference.
