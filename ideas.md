# Ideas

## Conformer Generation From Restraints

Generating full atom coordinates from restraints is much harder than hydrogen
placement, but still feasible if the goal is a single reasonable conformer
rather than a conformer ensemble.

Hydrogen generation is mostly local: once heavy-atom positions and local
geometry are known, each H can usually be placed from a small neighborhood.
Full conformer generation is globally coupled. Bond lengths and angles are
local, but torsions, chirality, ring closure, fused rings, macrocycles, and
nonbonded contacts interact across the whole graph.

Difficulty:

- H placement from restraints: low, mostly deterministic local construction.
- Heavy-atom xyz from restraints: medium-hard for one plausible conformer.
- Good conformer search / ensemble generation: hard.

Pragmatic Gemmi-native route:

1. Build an internal-coordinate tree where possible.
2. Seed rigid fragments first:
   - small rings,
   - aromatic systems,
   - planar groups,
   - peptide/nucleic templates when applicable.
3. Place the remaining graph by walking bonds using bond/angle/torsion
   restraints.
4. Enforce chirality during placement.
5. Run restrained Cartesian minimization to clean up closure errors and
   clashes.
6. For difficult rings/macrocycles, add a fallback embedding/minimization step.

Main hard parts:

- ring closure,
- mutually consistent torsions,
- chirality preservation during refinement,
- avoiding bad nonbonded contacts,
- disconnected or underdetermined restraint sets.

Practical conclusion:

- A first version that generates one deterministic idealized conformer for
  well-restrained small molecules looks achievable.
- Competing with RDKit/AceDRG-style conformer generation is a much larger
  project.


## Reproducing The COD -> Tables Part Of AceDRG

Reproducing the COD -> tables pipeline looks moderately hard, but much easier
than reproducing all of AceDRG.

This part is mainly an offline statistical extraction pipeline, not a full
runtime chemistry engine. That makes it tractable if the goal is useful
Gemmi-native tables rather than exact AceDRG reproduction.

Main hard parts:

- reproducing atom typing closely enough,
- matching fragmentation/environment rules,
- cleaning COD input robustly,
- choosing statistical policies for clustering, outlier rejection, sigma,
  torsion period, and priorities,
- making regeneration deterministic and reproducible.

Comparatively easier parts:

- parsing COD CIFs with Gemmi,
- extracting bonds/angles/torsions once graph + typing are defined,
- aggregating observations into tables,
- serializing Gemmi-native table files.

Difficulty by target:

- Useful Gemmi-native tables from COD: feasible.
- Faithful AceDRG-compatible reproduction of tables: hard.
- Scientifically reasonable replacement for AceDRG tables: feasible with
  sustained work.

Pragmatic implementation order:

1. Rebuild bond and angle tables first.
2. Define Gemmi-native atom/environment typing explicitly.
3. Add torsions later, starting with simpler non-ring torsions.
4. Keep peptide/nucleic/sugar special tables as separate curated layers.
5. Compare coverage and distributions against AceDRG rather than trying to
   match exact keys from the start.

Practical conclusion:

- A first useful version is a medium-size project.
- Exact AceDRG-style reproduction is a significantly larger project.

## Planar-Core Solver For `chemcomp_xyz`

The current XYZ generator handles tree-like fragments, small rings, simple
bridges, and terminal branches reasonably well, but it still fails badly on a
repeated motif: planar degree-3 junction atoms inside conjugated systems.

Typical failures seen in `ccd/gemmi/*`:

- `A5B`: `C8` and `N4` are shared junction atoms between planar fragments.
- `A6E`: `C6-N` bridge in a conjugated fragment.
- `A4W`: `C03-N04-C05` planar linker.
- `A0M`, `AK2`, `AVD`, `A2W`, `ABB`, `AW1`: amide/imide/arylamine-like
  trivalent `N`/`C` junctions.
- `A7Y`, `A6X`, `AWT`, `A55`, `AG7`, `A9G`, `ACS`: conjugated planar junction
  carbons inside fused or linked planar systems.

These failures are not mainly about hydrogen placement or small local cleanup.
They happen because a junction atom is placed from one local fragment only,
while its other planar neighbors belong to another fragment that is solved
later. Once the two sides drift apart, post-hoc plane projection or chirality
repair cannot recover the missing bond geometry.

### Goal

Add a dedicated planar-core embedding stage before generic atom-by-atom growth.
The planar-core stage should solve connected planar/conjugated heavy-atom
subgraphs as units in 2D, then hand the embedded core to the existing 3D
pipeline for branch growth, ring cleanup, and hydrogen placement.

### Scope Of The First Version

Handle only heavy-atom planar cores. Ignore hydrogens at this stage.

A planar core is a connected heavy-atom subgraph where atoms are linked by
bonds and are co-constrained by one or more `_chem_comp_plane_atom` groups.
The first version does not need to infer aromaticity or conjugation from bond
orders alone; it can rely primarily on explicit plane restraints plus bonded
connectivity.

### Why Single Plane Groups Were Not Enough

Trying to seed one `_chem_comp_plane_atom` group at a time did not fix `A5B`.
The key problem is overlap.

For example, a junction atom like `C8` can belong to multiple plane groups, and
those groups jointly define one larger planar core. Solving only one plane at a
time still leaves articulation atoms underconstrained at the interface between
planes.

So the unit of work has to be a connected planar-core component, not an
individual plane group.

### Data Model

Add a lightweight planar-core descriptor built from `ChemComp` restraints:

- `atoms`: indices of heavy atoms in the core.
- `bonds`: core-core bond list with target lengths.
- `angles`: in-core angles that help define local shape.
- `planes`: source plane-restraint group ids contributing to the core.
- `anchors`: atoms on the boundary that connect to non-core neighbors.
- `articulation_atoms`: atoms shared by multiple plane groups or linking major
  parts of the core.

Construction algorithm:

1. Start from each `_chem_comp_plane_atom` group.
2. Remove hydrogens.
3. Build a graph where heavy atoms are nodes and two atoms are linked if they
   are bonded and appear together in at least one plane group.
4. Merge overlapping plane groups by connected components in that graph.
5. Keep only components with at least 4 heavy atoms.

This gives a set of disjoint planar-core candidates.

### Embedding Strategy

Each planar core should be embedded in its own local 2D coordinate frame before
any generic growth pass touches its atoms.

First version, pragmatic route:

1. Pick a spanning tree of the core graph.
2. Choose one seed bond and place its two atoms on the x-axis using the target
   bond length.
3. Place a third atom using a target angle if available.
4. Grow the rest of the spanning tree in 2D using bond lengths and angles.
5. When an atom has multiple already placed core neighbors, choose its 2D
   position by satisfying multiple bond constraints, using angles only to break
   ties.
6. After all core atoms are placed, run a planar relaxation pass that averages
   inconsistent local placements and improves closure errors.

The result is a 2D embedding in a shared plane. It can then be lifted directly
into 3D with `z = 0` in the core-local frame.

### Solving Junction Atoms

The main reason to add this solver is to handle atoms such as:

- trivalent planar `N` with two carbon neighbors and one substituent,
- trivalent planar `C` shared between two conjugated fragments,
- atoms shared by two overlapping plane groups.

For such atoms, the local placement rule must prefer multi-neighbor solving:

- If 2 core neighbors are already placed: intersect two bond-distance
  constraints in the common plane.
- If 3 core neighbors are already placed: solve by least-squares position in
  the common plane, weighted by bond-length targets, then refine with angle
  targets if available.

This is different from the current one-center growth logic and is the core
reason the planar solver needs to exist as a separate stage.

### Relation To Existing `chemcomp_xyz` Pipeline

Recommended order:

1. Detect planar cores.
2. Embed planar cores as 2D heavy-atom fragments.
3. Embed non-planar ring seeds / other heavy-atom seeds.
4. Grow bridge atoms between already placed cores.
5. Grow remaining acyclic heavy branches.
6. Run current cleanup passes:
   - plane enforcement,
   - small-ring regularization,
   - safe chirality repair.
7. Add hydrogens last.

This keeps the current successful machinery for branches and hydrogens, while
moving the weak point (planar conjugated cores) into a dedicated stage.

### Acceptance Criteria

The first version does not need to solve all difficult ligands. It should be
judged on whether it removes the catastrophic failures for the dominant motif.

Useful target set for regression checks:

- `A5B`
- `A6E`
- `A4W`
- `A0M`
- `AK2`
- `AVD`
- `A2W`
- `ABB`
- `AW1`

Practical success criteria:

- no missing heavy-atom coordinates for these cases,
- no bond outliers larger than about `5 esd` for the worst previously broken
  planar-junction bonds,
- no regression in already improved cases such as `ABP` and `ALB`.

### Implementation Notes

- Keep it C++11-compatible.
- Avoid introducing another large pile of local repair heuristics into
  `chemcomp_xyz.cpp`; add a distinct helper stage instead.
- Do not try to solve full conformer generation here. This is only a planar
  heavy-core embedding stage.
- If the first version needs an iterative relaxation, keep it deterministic and
  lightweight.

### Suggested Incremental Steps

1. Add planar-core detection and debugging output only.
2. Add a 2D spanning-tree embedder for one connected planar core.
3. Add multi-neighbor placement for already placed core neighbors.
4. Run it only on a small whitelist in a script until it is stable.
5. Integrate it into `generate_chemcomp_xyz_from_restraints()`.
6. Expand regression checks to the broader `ccd/gemmi/*` corpus.
