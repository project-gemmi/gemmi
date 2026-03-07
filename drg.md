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

Remaining candidates:

1. Revisit output group labeling policy vs torsion mode policy.  
   Ensure `_chem_comp.group` and torsion-mode decisions remain chemically consistent.

2. Split explicit "compatibility mode" from default chemistry-first mode.  
   Keep strict AceDRG-emulation behavior behind a dedicated switch where still needed.
