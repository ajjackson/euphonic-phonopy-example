# Euphonic - phonopy interaction via Phonopy API

A demonstration script of Phonopy data loaded and converted to
Euphonic force constants without using the Euphonic phonopy file reader.

Usually Euphonic reads Phonopy output files directly without using the
Phonopy library, in order to avoid an unnecessary dependency. However,
for some workflows it makes sense to instantiate a Phonopy object
(e.g. to calculate the force constants) and convert to Euphonic
ForceConstants (e.g. for property calculations) in the same Python
session.

This repository is set up to use with `uv`, e.g.

```
uv run convert.py INPUT OUTPUT
```

will install dependencies to a clean virtual environment as required and run the format converter.
