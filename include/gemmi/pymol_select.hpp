// Copyright Global Phasing Ltd.
//
// PyMOL Selection syntax

#ifndef GEMMI_PYMO_SEL_HPP_
#define GEMMI_PYMO_SEL_HPP_

#include "flat.hpp"
#include "third_party/tao/pegtl.hpp" // IWYU pragma: keep

#include <string>
#include <vector>
#include <memory>
#include <iostream>

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
    std::string name;
    explicit ChainNode(std::string s) : name(std::move(s)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return a.chain_id == name;
    }
};

struct ResnNode : Node {
    std::string resn;
    explicit ResnNode(std::string s) : resn(std::move(s)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return a.residue_name == resn;
    }
};

struct AtomNameNode : Node {
    std::string name;
    explicit AtomNameNode(std::string s) : name(std::move(s)) {}
    bool match(const gemmi::FlatAtom& a) const override {
        return a.atom_name == name;
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

namespace p = tao::pegtl;

// --- State ---
struct State {
    std::vector<std::unique_ptr<psimpl::Node>> stack;
};

// --- Helpers ---
struct ws : p::star<p::space> {};
struct sep : p::plus<p::space> {}; // mandatory separator

// --- Values ---
struct integer : p::seq<p::opt<p::one<'-'>>, p::plus<p::digit>> {};
struct identifier : p::identifier {}; // standard C-style identifier
struct atom_name_str : p::plus<p::alnum> {}; // Allow 1+2 for names? keeping simple for now

// --- Keywords ---
// Using istring for case-insensitive matching
struct kw_chain : p::istring<'c','h','a','i','n'> {};
struct kw_resn  : p::istring<'r','e','s','n'> {};
struct kw_resi  : p::istring<'r','e','s','i'> {};
struct kw_name  : p::istring<'n','a','m','e'> {};
struct kw_alt   : p::istring<'a','l','t'> {};

struct kw_and   : p::istring<'a','n','d'> {};
struct kw_or    : p::istring<'o','r'> {};
struct kw_not   : p::istring<'n','o','t'> {};

// --- Property Rules ---

// Chain: chain A
struct val_chain : identifier {};
struct rule_chain : p::seq<kw_chain, sep, val_chain> {};

// Resn: resn ALA
struct val_resn : atom_name_str {};
struct rule_resn : p::seq<kw_resn, sep, val_resn> {};

// Name: name CA
struct val_name : atom_name_str {};
struct rule_name : p::seq<kw_name, sep, val_name> {};

// Alt: alt A
struct val_alt : p::alnum {}; // single char
struct rule_alt : p::seq<kw_alt, sep, val_alt> {};

// Resi: resi 100 OR resi 100-200
struct val_resi_range : p::seq<integer, p::one<'-'>, integer> {};
struct val_resi_single : integer {};
struct rule_resi : p::seq<kw_resi, sep, p::sor<val_resi_range, val_resi_single>> {};

// Combined Property
struct property : p::sor<
    rule_chain,
    rule_resn,
    rule_resi,
    rule_name,
    rule_alt
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

// --- Leaf Actions ---

template<> struct action<val_chain> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.stack.push_back(std::make_unique<psimpl::ChainNode>(in.string()));
    }
};

template<> struct action<val_resn> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.stack.push_back(std::make_unique<psimpl::ResnNode>(in.string()));
    }
};

template<> struct action<val_name> {
    template<typename Input> static void apply(const Input& in, State& s) {
        s.stack.push_back(std::make_unique<psimpl::AtomNameNode>(in.string()));
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
        // Input is "100-200". We need to split it.
        std::string str = in.string();
        // Handle negative numbers logic if needed, but keep simple for "int-int"
        // If the first int is negative, find might hit index 0.
        // A robust split for "Start-End" requires parsing context or re-parsing.

        // Quick hack for this action: find the dash that separates the two numbers.
        // Since `integer` rule consumes greedy digits, the PEGTL match guarantees structure.
        // But standard string split is messy with potential negatives "-5--2".
        // Let's assume positive ranges for the "simple" single file implementation
        // OR re-parse using simple string logic.

        size_t split_pos = str.find('-', 1); // Skip potential leading negative sign
        int v1 = std::stoi(str.substr(0, split_pos));
        int v2 = std::stoi(str.substr(split_pos + 1));
        s.stack.push_back(std::make_unique<psimpl::ResiRangeNode>(v1, v2));
    }
};

// --- Logic Actions ---

// NOT
template<> struct action<rule_not> {
    static void apply0(State& s) {
        // The child was just pushed by 'factor'.
        auto child = std::move(s.stack.back());
        s.stack.pop_back();
        s.stack.push_back(std::make_unique<psimpl::NotNode>(std::move(child)));
    }
};

// AND
template<> struct action<and_rest> {
    static void apply0(State& s) {
        // Stack: [LHS, RHS] -> [AndNode(LHS, RHS)]
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

// Applies the selection to a Gemmi structure
inline std::vector<gemmi::FlatAtom*> select_atoms(gemmi::Structure& st, const std::string& query) {
    auto root = compile_pymol_selection(query);
    gemmi::FlatStructure fs(st);
    std::vector<gemmi::FlatAtom*> result;

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
