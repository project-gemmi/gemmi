# pattern file for nanobind stubgen.py

# Here, we intentionally use narrowed arguments to prevent some mistakes.

FTransform.apply:
    def apply(self, arg: Fractional, /) -> Fractional: ...  # type: ignore[override]

Fractional.__(add|sub)__:
    def __\1__(self, arg: Fractional, /) -> Fractional: ...  # type: ignore[override]

Position.__(add|sub|iadd|isub)__:
    def __\1__(self, arg: Position, /) -> Position: ...  # type: ignore[override]

# Here, ResidueGroup group inherits implementation from ResidueSpan,
# but the items are accessed differently.

ResidueGroup.__getitem__:
    @overload  # type: ignore[override]
    def __getitem__(self, index: int) -> Residue: ...

    @overload
    def __getitem__(self, name: str) -> Residue: ...

ResidueGroup.__delitem__:
    def __delitem__(self, name: str) -> None: ...  # type: ignore[override]

# atm, IntFlag-like enums generate wrong __str__, disable it for now
ChemCompModel.__(str|repr)__:
    # __str__ and __repr__ disabled in stubs.pat
