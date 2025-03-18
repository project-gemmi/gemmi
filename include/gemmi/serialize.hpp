// Copyright Global Phasing Ltd.
//
// Binary serialization for Structure (as well as Model, UnitCell, etc).
//
// Based on zpp::serializer, include third_party/serializer.h first.

#ifndef GEMMI_SERIALIZE_HPP_
#define GEMMI_SERIALIZE_HPP_

#include "model.hpp"
#include "cifdoc.hpp"

#define SERIALIZE(Struct, ...) \
template <typename Archive> \
void serialize(Archive& archive, Struct& o) { archive(__VA_ARGS__); } \
template <typename Archive> \
void serialize(Archive& archive, const Struct& o) { archive(__VA_ARGS__); }

#define SERIALIZE_P(Struct, Parent, ...) \
template <typename Archive> \
void serialize(Archive& archive, Struct& o) \
 { archive(static_cast<Parent&>(o), __VA_ARGS__); } \
template <typename Archive> \
void serialize(Archive& archive, const Struct& o) \
 { archive(static_cast<const Parent&>(o), __VA_ARGS__); }

#define SERIALIZE_T1(Struct, Typename, ...) \
template <typename Archive, Typename T> \
void serialize(Archive& archive, Struct<T>& o) { archive(__VA_ARGS__); } \
template <typename Archive, Typename T> \
void serialize(Archive& archive, const Struct<T>& o) { archive(__VA_ARGS__); }

namespace gemmi {

SERIALIZE_T1(OptionalInt, int, o.value)

//SERIALIZE(Element, o.elem) is ambiguous because of El->Element conversion
template <typename Archive>
void serialize(Archive& archive, Element& o) { archive((unsigned char&)o.elem); }
template <typename Archive>
void serialize(Archive& archive, const Element& o) { archive((unsigned char)o.elem); }

SERIALIZE_T1(Vec3_, typename, o.x, o.y, o.z)

SERIALIZE_T1(SMat33, typename, o.u11, o.u22, o.u33, o.u12, o.u13, o.u23)

SERIALIZE(Mat33, o.a)

SERIALIZE(Transform, o.mat, o.vec)

SERIALIZE(NcsOp, o.id, o.given, o.tr)

SERIALIZE(UnitCellParameters, o.a, o.b, o.c, o.alpha, o.beta, o.gamma)

SERIALIZE_P(UnitCell, UnitCellParameters,
          o.orth, o.frac, o.volume,
          o.ar, o.br, o.cr, o.cos_alphar, o.cos_betar, o.cos_gammar,
          o.explicit_matrices, o.cs_count, o.images)

SERIALIZE(SeqId, o.num, o.icode)

SERIALIZE(AtomAddress, o.chain_name, o.res_id, o.atom_name, o.altloc)

SERIALIZE(Metadata, o.authors, o.experiments, o.crystals, o.refinement,
          o.software, o.solved_by, o.starting_model, o.remark_300_detail)

SERIALIZE(SoftwareItem, o.name, o.version, o.date, o.description,
          o.contact_author, o.contact_author_email, o.classification)

SERIALIZE(ReflectionsInfo, o.resolution_high, o.resolution_low, o.completeness,
          o.redundancy, o.r_merge, o.r_sym, o.mean_I_over_sigma)

SERIALIZE(ExperimentInfo, o.method, o.number_of_crystals, o.unique_reflections,
          o.reflections, o.b_wilson, o.shells, o.diffraction_ids)

SERIALIZE(DiffractionInfo, o.id, o.temperature, o.source, o.source_type,
          o.synchrotron, o.beamline, o.wavelengths, o.scattering_type,
          o.mono_or_laue, o.monochromator, o.collection_date, o.optics,
          o.detector, o.detector_make)

SERIALIZE(CrystalInfo, o.id, o.description, o.ph, o.ph_range, o.diffractions)

SERIALIZE(TlsGroup, o.num_id, o.id, o.selections, o.origin, o.T, o.L, o.S)

SERIALIZE(TlsGroup::Selection, o.chain, o.res_begin, o.res_end, o.details)

SERIALIZE(BasicRefinementInfo, o.resolution_high, o.resolution_low,
          o.completeness, o.reflection_count, o.work_set_count, o.rfree_set_count,
          o.r_all, o.r_work, o.r_free,
          o.cc_fo_fc_work, o.cc_fo_fc_free, o.fsc_work, o.fsc_free,
          o.cc_intensity_work, o.cc_intensity_free)

SERIALIZE_P(RefinementInfo, BasicRefinementInfo, o.id,
            o.cross_validation_method, o.rfree_selection_method,
            o.bin_count, o.bins, o.mean_b, o.aniso_b,
            o.luzzati_error, o.dpi_blow_r, o.dpi_blow_rfree,
            o.dpi_cruickshank_r, o.dpi_cruickshank_rfree,
            o.restr_stats, o.tls_groups, o.remarks)

SERIALIZE(RefinementInfo::Restr, o.name, o.count, o.weight, o.function, o.dev_ideal)

SERIALIZE(Entity, o.name, o.subchains, o.entity_type, o.polymer_type,
          o.reflects_microhetero, o.dbrefs, o.sifts_unp_acc, o.full_sequence)

SERIALIZE(Entity::DbRef, o.db_name, o.accession_code, o.id_code,
          o.isoform, o.seq_begin, o.seq_end, o.db_begin, o.db_end,
          o.label_seq_begin, o.label_seq_end)

SERIALIZE(SiftsUnpResidue, o.res, o.acc_index, o.num)

SERIALIZE(Connection, o.name, o.link_id, o.type, o.asu,
          o.partner1, o.partner2, o.reported_distance, o.reported_sym)

SERIALIZE(CisPep, o.partner_c, o.partner_n, o.model_num, o.only_altloc, o.reported_angle)

SERIALIZE(ModRes, o.chain_name, o.res_id, o.parent_comp_id, o.mod_id, o.details)

SERIALIZE(Helix, o.start, o.end, o.pdb_helix_class, o.length)

SERIALIZE(Sheet, o.name, o.strands)

SERIALIZE(Sheet::Strand, o.start, o.end, o.hbond_atom2, o.hbond_atom1, o.sense, o.name)

SERIALIZE(Assembly, o.name, o.author_determined, o.software_determined,
          o.special_kind, o.oligomeric_count, o.oligomeric_details,
          o.software_name, o.absa, o.ssa, o.more, o.generators)

SERIALIZE(Assembly::Operator, o.name, o.type, o.transform)

SERIALIZE(Assembly::Gen, o.chains, o.subchains, o.operators)

SERIALIZE(ResidueId, o.seqid, o.segment, o.name)

SERIALIZE(Atom, o.name, o.altloc, o.charge, o.element, o.calc_flag,
          o.flag, o.tls_group_id, o.serial, o.fraction, o.pos,
          o.occ, o.b_iso, o.aniso)

SERIALIZE_P(Residue, ResidueId,
            o.subchain, o.entity_id, o.label_seq, o.entity_type,
            o.het_flag, o.flag, o.sifts_unp, o.group_idx, o.atoms)

SERIALIZE(Chain, o.name, o.residues)

SERIALIZE(Model, o.num, o.chains)

SERIALIZE(Structure, o.name, o.cell, o.spacegroup_hm, o.models,
          o.ncs, o.entities, o.connections, o.cispeps, o.mod_residues,
          o.helices, o.sheets, o.assemblies, o.conect_map, o.meta,
          o.input_format, o.has_d_fraction, o.ter_status,
          o.has_origx, o.origx, o.info, o.shortened_ccd_codes,
          o.raw_remarks, o.resolution)


namespace cif {

SERIALIZE(Loop, o.tags, o.values)
SERIALIZE(Block, o.name, o.items)
SERIALIZE(Document, o.source, o.blocks)

template <typename Archive>
void serialize(Archive& archive, Item& o) {
  archive(o.type, o.line_number);
  switch (o.type) {
    case ItemType::Pair:
    case ItemType::Comment: new(&o.pair) cif::Pair; archive(o.pair); break;
    case ItemType::Loop: new(&o.loop) cif::Loop; archive(o.loop); break;
    case ItemType::Frame: new(&o.frame) cif::Block; archive(o.frame); break;
    case ItemType::Erased: break;
  }
}
template <typename Archive>
void serialize(Archive& archive, const Item& o) {
  archive(o.type, o.line_number);
  switch (o.type) {
    case ItemType::Pair:
    case ItemType::Comment: archive(o.pair); break;
    case ItemType::Loop: archive(o.loop); break;
    case ItemType::Frame: archive(o.frame); break;
    case ItemType::Erased: break;
  }
}

} // namespace cif
} // namespace gemmi

#undef SERIALIZE
#undef SERIALIZE_P
#undef SERIALIZE_T1

#endif
