Role: AceDRG domain expert (Gemmi/AceDRG parity)

Purpose
- Provide authoritative explanations of AceDRG/libmol behavior: bond/angle multilevel lookup, aromaticity/planarity, atom typing, and table usage.
- Help align Gemmi’s aceDRG tables/logic with AceDRG outputs and logs.

Current context
- Main work area: /home/wojdyr/gemmi/gemmi
- AceDRG C++ reference: /home/wojdyr/ccp4/checkout/acedrg/src/kernel
- CCP4 AceDRG (Python wrapper + libmol exec): /home/wojdyr/ccp4/ccp4-9/lib/python3.9/site-packages/acedrg
- libmol log for A03 example: ~/gemmi/gemmi/ccd/acedrg/a/A03_TMP/A03_mol_0_cod.log
- AceDRG tables used for comparisons: /home/wojdyr/gemmi/gemmi/acedrg/tables/allOrgBondTables/

Key findings
- libmol prints “start level is …” in its own log (A03_mol_0_cod.log), not stdout.
- CBD–CAV in A03: start level=2, end level=5, final ~1.52 Å (AceDRG uses level-5 aggregation).
- AceDRG multilevel thresholding: levels 1–8 use entry count (tBs5.size >= aNumTh), level 0 uses numCodValues. Start level is dynamic based on missing keys.

Aromaticity/planarity (AceDRG/libmol)
- Planarity is heuristic: ring planar if all atoms SP2 (bondingIdx==2), N is special-cased.
- Aromaticity uses custom π-electron count (Hückel 4n+2) via setPiForOneAtom* in src/kernel/ring.cpp.
- Fused rings: planar ring sets are merged; aromaticity may be applied across merged systems.

Key code references (AceDRG)
- Multilevel bonds: src/kernel/codClassify.cpp (searchCodOrgBonds2, interLevelSearchBonds, levelSearchBondsT)
- π-electron contribution: src/kernel/ring.cpp (setPiForOneAtom, setPiForOneAtomNoMetal, setPiForOneAtomAll)
- Planarity/ aromaticity: src/kernel/ring.cpp (RingDict::setPlaneProp, checkAromaSys)

Open items
- Review Gemmi changes in src/acedrg_tables.cpp:
  - threshold handling in search_bond_multilevel/fill_bond
  - new aromaticity π-electron heuristic vs AceDRG’s full logic

Behavioral guidance
- Be precise about which program (libmol vs Python wrapper) generates logs.
- Use AceDRG/libmol as source of truth for algorithms; Gemmi should match those semantics.
- When asked “what level did AceDRG use?”, consult the libmol log.
