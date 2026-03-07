# AceDRG Compatibility Cleanup Status

Items completed in Gemmi:

1. Removed descriptor-presence gating.  
   `_pdbx_chem_comp_descriptor` no longer controls peptide/nucleic torsion correction flow.

2. Enabled nucleic torsion-table replacement for nucleic components.  
   DNA/RNA-like types now apply nucleic overrides directly.

3. Removed strict peptide backbone naming/topology gate.  
   Peptide override use is no longer blocked by exact N/CA/C/O/OXT and H-name checks.

4. Replaced branchy AceDRG-like torsion candidate selection with scored selection.  
   Selection now prefers chemically informative termini (non-H/ring-aware) and rejects self-torsions.

5. Made pyranose chair enforcement conditional on coordinate completeness.  
   With complete coordinates, coordinate-driven torsions are kept.

6. Added explicit compatibility mode switch.  
   `GEMMI_ACE_COMPAT=1` now enables the AceDRG-like peptide gate, torsion candidate
   selector, and pyranose-chair behavior; default mode stays chemistry-first.

7. Harmonized output group inference with type metadata.  
   When `cc.group` is missing, infer it from `type_or_group` (and residue-class fallback)
   to keep group labeling consistent with torsion mode intent.

Remaining candidates:

1. Review whether any other `*_like_acedrg` branches should be compat-only.  
   Keep default behavior chemistry-driven unless strict emulation is requested.
