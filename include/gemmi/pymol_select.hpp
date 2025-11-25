// Copyright Global Phasing Ltd.
//
// PyMOL Selection syntax

#ifndef GEMMI_PYMO_SEL_HPP_
#define GEMMI_PYMO_SEL_HPP_

#include "flat.hpp"
#include "glob.hpp"  // for glob_match
#include "third_party/tao/pegtl.hpp" // IWYU pragma: keep

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <cstring>

// ============================================================================
// PHASE 1: Abstract Syntax Tree (AST)
// ============================================================================


namespace gemmi {

namespace psimpl {

struct Node {
    virtual ~Node() = default;
    virtual bool match(const gemmi::FlatAtom& a) const = 0;
};

// --- Logic Nodes ---

struct AndNode : Node {
    std::unique_ptr<Node> left, right;
    bool match(const gemmi::FlatAtom& a) const override {
        return left->match(a) && right->match(a);
    }
};

struct OrNode : Node {
    std::unique_ptr<Node> left, right;
    bool match(const gemmi::FlatAtom& a) const override {
        return left->match(a) || right->match(a);
    }
};

struct NotNode : Node {
    std::unique_ptr<Node> child;
    explicit NotNode(std::unique_ptr<Node> c) : child(std::move(c)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return !child->match(a);
    }
};

// --- Property Nodes ---

struct ChainNode : Node {
    std::vector<std::string> names;
    explicit ChainNode(std::vector<std::string> v) : names(std::move(v)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        for (const auto& n : names)
            if (glob_match(n, a.chain_id)) return true;
        return false;
    }
};

struct ResnNode : Node {
    std::vector<std::string> names;
    explicit ResnNode(std::vector<std::string> v) : names(std::move(v)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        for (const auto& n : names)
            if (glob_match(n, a.residue_name)) return true;
        return false;
    }
};

struct AtomNameNode : Node {
    std::vector<std::string> names;
    explicit AtomNameNode(std::vector<std::string> v) : names(std::move(v)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        for (const auto& n : names)
            if (glob_match(n, a.atom_name)) return true;
        return false;
    }
};

struct AltLocNode : Node {
    char alt;
    explicit AltLocNode(char c) : alt(c) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return a.altloc == alt;
    }
};

struct ResiRangeNode : Node {
    int min, max;
    ResiRangeNode(int a, int b) : min(a), max(b) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return *a.seq_id.num >= min && *a.seq_id.num <= max;
    }
};

struct ElementNode : Node {
    std::vector<Element> elems;
    explicit ElementNode(std::vector<Element> v) : elems(std::move(v)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        for (const auto& e : elems)
            if (a.element == e) return true;
        return false;
    }
};

struct HetatmNode : Node {
    bool hetatm;  // true = hetatm, false = not hetatm (i.e., ATOM)
    explicit HetatmNode(bool h) : hetatm(h) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return hetatm ? (a.het_flag == 'H') : (a.het_flag == 'A');
    }
};

struct EntityTypeNode : Node {
    EntityType etype;
    explicit EntityTypeNode(EntityType e) : etype(e) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return a.entity_type == etype;
    }
};

struct HydrogenNode : Node {
    bool match(const gemmi::FlatAtom& a) const override {
        return a.element == El::H || a.element == El::D;
    }
};

enum class CompareOp { LT, LE, GT, GE, EQ, NE };

struct BfactorNode : Node {
    CompareOp op;
    float value;
    BfactorNode(CompareOp o, float v) : op(o), value(v) {}
    bool match(const gemmi::FlatAtom& a) const override {
        switch (op) {
            case CompareOp::LT: return a.b_iso < value;
            case CompareOp::LE: return a.b_iso <= value;
            case CompareOp::GT: return a.b_iso > value;
            case CompareOp::GE: return a.b_iso >= value;
            case CompareOp::EQ: return a.b_iso == value;
            case CompareOp::NE: return a.b_iso != value;
        }
        return false;
    }
};

struct OccupancyNode : Node {
    CompareOp op;
    float value;
    OccupancyNode(CompareOp o, float v) : op(o), value(v) {}
    bool match(const gemmi::FlatAtom& a) const override {
        switch (op) {
            case CompareOp::LT: return a.occ < value;
            case CompareOp::LE: return a.occ <= value;
            case CompareOp::GT: return a.occ > value;
            case CompareOp::GE: return a.occ >= value;
            case CompareOp::EQ: return a.occ == value;
            case CompareOp::NE: return a.occ != value;
        }
        return false;
    }
};

struct BackboneNode : Node {
    bool match(const gemmi::FlatAtom& a) const override {
        // Standard protein backbone atoms
        return std::strcmp(a.atom_name, "CA") == 0 ||
               std::strcmp(a.atom_name, "C") == 0 ||
               std::strcmp(a.atom_name, "N") == 0 ||
               std::strcmp(a.atom_name, "O") == 0;
    }
};

struct SidechainNode : Node {
    bool match(const gemmi::FlatAtom& a) const override {
        // Sidechain = not backbone and not hydrogen
        return std::strcmp(a.atom_name, "CA") != 0 &&
               std::strcmp(a.atom_name, "C") != 0 &&
               std::strcmp(a.atom_name, "N") != 0 &&
               std::strcmp(a.atom_name, "O") != 0 &&
               a.element != El::H && a.element != El::D;
    }
};

struct AllNode : Node {
    bool match(const gemmi::FlatAtom&) const override { return true; }
};

namespace p = tao::pegtl;

// --- State ---
struct State {
    std::vector<std::unique_ptr<psimpl::Node>> stack;
    std::vector<std::string> string_list;  // temp storage for building value lists
    CompareOp current_op = CompareOp::EQ;
};

// --- Helpers ---
struct ws : p::star<p::space> {};
struct sep : p::plus<p::space> {}; // mandatory separator

// --- Values ---
struct integer : p::seq<p::opt<p::one<'-'>>, p::plus<p::digit>> {};
struct float_num : p::seq<p::opt<p::one<'-'>>, p::plus<p::digit>,
                          p::opt<p::seq<p::one<'.'>, p::star<p::digit>>>> {};
// Allow wildcards in identifiers and names
struct wildcard_char : p::one<'*', '?'> {};
struct identifier : p::plus<p::sor<p::alnum, p::one<'_'>, wildcard_char>> {};
struct atom_name_str : p::plus<p::sor<p::alnum, wildcard_char>> {};
struct element_str : p::seq<p::upper, p::opt<p::lower>> {}; // e.g., C, Ca, Fe (no wildcards for elements)

// --- Comparison operators ---
struct op_le : p::string<'<','='> {};
struct op_ge : p::string<'>','='> {};
struct op_ne : p::sor<p::string<'!','='>, p::string<'<','>'>> {};
struct op_lt : p::one<'<'> {};
struct op_gt : p::one<'>'> {};
struct op_eq : p::one<'='> {};
struct compare_op : p::sor<op_le, op_ge, op_ne, op_lt, op_gt, op_eq> {};

// --- Keywords ---
// Using istring for case-insensitive matching
struct kw_chain : p::istring<'c','h','a','i','n'> {};
struct kw_resn  : p::istring<'r','e','s','n'> {};
struct kw_resi  : p::istring<'r','e','s','i'> {};
struct kw_name  : p::istring<'n','a','m','e'> {};
struct kw_alt   : p::istring<'a','l','t'> {};
struct kw_elem  : p::istring<'e','l','e','m'> {};
struct kw_b     : p::istring<'b'> {};
struct kw_q     : p::istring<'q'> {};

struct kw_and   : p::istring<'a','n','d'> {};
struct kw_or    : p::istring<'o','r'> {};
struct kw_not   : p::istring<'n','o','t'> {};

// Stand-alone keywords (no arguments)
struct kw_hetatm    : p::istring<'h','e','t','a','t','m'> {};
struct kw_polymer   : p::istring<'p','o','l','y','m','e','r'> {};
struct kw_solvent   : p::istring<'s','o','l','v','e','n','t'> {};
struct kw_water     : p::istring<'w','a','t','e','r'> {};
struct kw_hydrogens : p::istring<'h','y','d','r','o','g','e','n','s'> {};
struct kw_h_dot     : p::istring<'h','.'> {};
struct kw_backbone  : p::istring<'b','a','c','k','b','o','n','e'> {};
struct kw_sidechain : p::istring<'s','i','d','e','c','h','a','i','n'> {};
struct kw_all       : p::istring<'a','l','l'> {};

// --- Property Rules ---

// Chain: chain A or chain A+B+C
struct val_chain_item : identifier {};
struct val_chain_list : p::list<val_chain_item, p::one<'+'>> {};
struct rule_chain : p::seq<kw_chain, sep, val_chain_list> {};

// Resn: resn ALA or resn ALA+GLY+VAL
struct val_resn_item : atom_name_str {};
struct val_resn_list : p::list<val_resn_item, p::one<'+'>> {};
struct rule_resn : p::seq<kw_resn, sep, val_resn_list> {};

// Name: name CA or name CA+CB+N
struct val_name_item : atom_name_str {};
struct val_name_list : p::list<val_name_item, p::one<'+'>> {};
struct rule_name : p::seq<kw_name, sep, val_name_list> {};

// Alt: alt A
struct val_alt : p::alnum {}; // single char
struct rule_alt : p::seq<kw_alt, sep, val_alt> {};

// Resi: resi 100 OR resi 100-200
struct val_resi_range : p::seq<integer, p::one<'-'>, integer> {};
struct val_resi_single : integer {};
struct rule_resi : p::seq<kw_resi, sep, p::sor<val_resi_range, val_resi_single>> {};

// Elem: elem C or elem C+N+O
struct val_elem_item : element_str {};
struct val_elem_list : p::list<val_elem_item, p::one<'+'>> {};
struct rule_elem : p::seq<kw_elem, sep, val_elem_list> {};

// B-factor: b > 50, b < 20, b = 0
struct val_b_compare : float_num {};
struct rule_b : p::seq<kw_b, ws, compare_op, ws, val_b_compare> {};

// Occupancy: q < 1, q > 0.5
struct val_q_compare : float_num {};
struct rule_q : p::seq<kw_q, ws, compare_op, ws, val_q_compare> {};

// Stand-alone keywords
struct rule_hetatm    : kw_hetatm {};
struct rule_polymer   : kw_polymer {};
struct rule_solvent   : kw_solvent {};
struct rule_water     : kw_water {};
struct rule_hydrogens : p::sor<kw_hydrogens, kw_h_dot> {};
struct rule_backbone  : kw_backbone {};
struct rule_sidechain : kw_sidechain {};
struct rule_all       : kw_all {};

// Combined Property
struct property : p::sor<
    rule_chain,
    rule_resn,
    rule_resi,
    rule_name,
    rule_alt,
    rule_elem,
    rule_b,
    rule_q,
    rule_hetatm,
    rule_polymer,
    rule_solvent,
    rule_water,
    rule_hydrogens,
    rule_backbone,
    rule_sidechain,
    rule_all
> {};

// --- Boolean Logic Rules ---

struct expression; // forward decl

struct parens : p::seq<p::one<'('>, ws, expression, ws, p::one<')'>> {};

// Factor: NOT factor | parens | property
struct not_factor : p::seq<kw_not, sep, p::seq<expression>> {}; // simplified recursion
// Actually, to handle precedence properly with PEGTL without left-recursion:
// factor = (NOT ws factor) | parens | property
struct factor;
struct rule_not : p::seq<kw_not, ws, factor> {};
struct factor : p::sor<rule_not, parens, property> {};

// Term (AND)
// To allow easy "reduce" actions, we explicitly name the sequence
struct and_rest : p::seq<ws, kw_and, ws, factor> {};
struct term : p::seq<factor, p::star<and_rest>> {};

// Expression (OR)
struct or_rest : p::seq<ws, kw_or, ws, term> {};
struct expression : p::seq<term, p::star<or_rest>> {};

// Root
struct grammar : p::must<ws, expression, ws, p::eof> {};

// ============================================================================
// PHASE 3: Actions
// ============================================================================

template<typename Rule>
struct action : p::nothing<Rule> {};

// --- List item actions (accumulate into string_list) ---

template<> struct action<val_chain_item> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.string_list.push_back(in.string());
    }
};

template<> struct action<val_resn_item> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.string_list.push_back(in.string());
    }
};

template<> struct action<val_name_item> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.string_list.push_back(in.string());
    }
};

template<> struct action<val_elem_item> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.string_list.push_back(in.string());
    }
};

// --- Rule actions (create nodes from accumulated lists) ---

template<> struct action<rule_chain> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::ChainNode>(std::move(s.string_list)));
        s.string_list.clear();
    }
};

template<> struct action<rule_resn> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::ResnNode>(std::move(s.string_list)));
        s.string_list.clear();
    }
};

template<> struct action<rule_name> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::AtomNameNode>(std::move(s.string_list)));
        s.string_list.clear();
    }
};

template<> struct action<rule_elem> {
    static void apply0(State& s) {
        std::vector<Element> elems;
        for (const auto& str : s.string_list)
            elems.push_back(Element(str));
        s.stack.push_back(std::make_unique<psimpl::ElementNode>(std::move(elems)));
        s.string_list.clear();
    }
};

template<> struct action<val_alt> {
    template<typename Input> static void apply(const Input& in, State& s) {
        std::string str = in.string();
        char c = str.empty() ? ' ' : str[0];
        s.stack.push_back(std::make_unique<psimpl::AltLocNode>(c));
    }
};

template<> struct action<val_resi_single> {
    template<typename Input> static void apply(const Input& in, State& s) {
        int val = std::stoi(in.string());
        s.stack.push_back(std::make_unique<psimpl::ResiRangeNode>(val, val));
    }
};

template<> struct action<val_resi_range> {
    template<typename Input> static void apply(const Input& in, State& s) {
        std::string str = in.string();
        size_t split_pos = str.find('-', 1); // Skip potential leading negative sign
        int v1 = std::stoi(str.substr(0, split_pos));
        int v2 = std::stoi(str.substr(split_pos + 1));
        s.stack.push_back(std::make_unique<psimpl::ResiRangeNode>(v1, v2));
    }
};

// --- Comparison operator actions ---

template<> struct action<op_lt> {
    static void apply0(State& s) { s.current_op = CompareOp::LT; }
};
template<> struct action<op_le> {
    static void apply0(State& s) { s.current_op = CompareOp::LE; }
};
template<> struct action<op_gt> {
    static void apply0(State& s) { s.current_op = CompareOp::GT; }
};
template<> struct action<op_ge> {
    static void apply0(State& s) { s.current_op = CompareOp::GE; }
};
template<> struct action<op_eq> {
    static void apply0(State& s) { s.current_op = CompareOp::EQ; }
};
template<> struct action<op_ne> {
    static void apply0(State& s) { s.current_op = CompareOp::NE; }
};

template<> struct action<val_b_compare> {
    template<typename Input> static void apply(const Input& in, State& s) {
        float val = std::stof(in.string());
        s.stack.push_back(std::make_unique<psimpl::BfactorNode>(s.current_op, val));
    }
};

template<> struct action<val_q_compare> {
    template<typename Input> static void apply(const Input& in, State& s) {
        float val = std::stof(in.string());
        s.stack.push_back(std::make_unique<psimpl::OccupancyNode>(s.current_op, val));
    }
};

// --- Stand-alone keyword actions ---

template<> struct action<rule_hetatm> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::HetatmNode>(true));
    }
};

template<> struct action<rule_polymer> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::EntityTypeNode>(EntityType::Polymer));
    }
};

template<> struct action<rule_solvent> {
    static void apply0(State& s) {
        // Solvent includes water
        s.stack.push_back(std::make_unique<psimpl::EntityTypeNode>(EntityType::Water));
    }
};

template<> struct action<rule_water> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::EntityTypeNode>(EntityType::Water));
    }
};

template<> struct action<rule_hydrogens> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::HydrogenNode>());
    }
};

template<> struct action<rule_backbone> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::BackboneNode>());
    }
};

template<> struct action<rule_sidechain> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::SidechainNode>());
    }
};

template<> struct action<rule_all> {
    static void apply0(State& s) {
        s.stack.push_back(std::make_unique<psimpl::AllNode>());
    }
};

// --- Logic Actions ---

// NOT
template<> struct action<rule_not> {
    static void apply0(State& s) {
        auto child = std::move(s.stack.back());
        s.stack.pop_back();
        s.stack.push_back(std::make_unique<psimpl::NotNode>(std::move(child)));
    }
};

// AND
template<> struct action<and_rest> {
    static void apply0(State& s) {
        auto rhs = std::move(s.stack.back()); s.stack.pop_back();
        auto lhs = std::move(s.stack.back()); s.stack.pop_back();

        auto node = std::make_unique<psimpl::AndNode>();
        node->left = std::move(lhs);
        node->right = std::move(rhs);
        s.stack.push_back(std::move(node));
    }
};

// OR
template<> struct action<or_rest> {
    static void apply0(State& s) {
        auto rhs = std::move(s.stack.back()); s.stack.pop_back();
        auto lhs = std::move(s.stack.back()); s.stack.pop_back();

        auto node = std::make_unique<psimpl::OrNode>();
        node->left = std::move(lhs);
        node->right = std::move(rhs);
        s.stack.push_back(std::move(node));
    }
};

} // namespace psimpl


// ============================================================================
// Public API
// ============================================================================

// Returns a compiled selection tree
inline std::unique_ptr<psimpl::Node> compile_pymol_selection(const std::string& selector) {
    psimpl::State state;
    tao::pegtl::memory_input<> in(selector, "");
    try {
        tao::pegtl::parse<psimpl::grammar, psimpl::action>(in, state);
        if (state.stack.empty()) return nullptr;
        return std::move(state.stack.back());
    }
    catch (const tao::pegtl::parse_error& e) {
        std::cerr << "Selection Parse Error: " << e.what() << std::endl;
        return nullptr;
    }
}

inline std::vector<const gemmi::FlatAtom*>
select_atoms(const gemmi::FlatStructure& fs, const std::string& query) {
    auto root = compile_pymol_selection(query);
    std::vector<const gemmi::FlatAtom*> result;

    if (root)
      for (auto& atom : fs.table) {
        if (root->match(atom)) {
          result.push_back(&atom);
        }
    }
    return result;
}

} // namespace gemmi

#endif
