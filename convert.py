from argparse import ArgumentParser
from pathlib import Path
import re

from euphonic import Crystal, ForceConstants, ureg
from euphonic.util import convert_fc_phases
import numpy as np
import phonopy


def get_parser() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument(
        "phonopy_yaml",
        type=Path,
        nargs="?",
        default=Path(__file__).parent / "data/Si-phonopy.yaml",
    )
    parser.add_argument(
        "output_json",
        type=Path,
        nargs="?",
        default=Path("euphonic-fc.json"),
    )
    return parser


def get_phonopy(phonopy_yaml: Path) -> phonopy.Phonopy:
    """Read phonopy.yaml with HDF5 force constants if available"""

    re_match = re.match(r"(?P<seedname>.+)-phonopy\.ya?ml", phonopy_yaml.name)
    if re_match:
        seedname = re_match.group("seedname")
        fc_file = phonopy_yaml.parent / f"{seedname}-force_constants.hdf5"
        if not fc_file.exists():
            fc_file = None
    else:
        fc_file = None

    phonon = phonopy.load(
        str(phonopy_yaml),
        force_constants_filename=str(fc_file),
    )
    return phonon


def phonopy_to_euphonic(phonon: phonopy.Phonopy) -> ForceConstants:
    """Convert from Phonopy to Euphonic force constant objects"""

    if phonon.nac_params:
        raise NotImplementedError(
            "Born/dielectric loading not yet implemented here; shouldn't be too hard!"
        )

    primitive = phonon.primitive
    supercell = phonon.supercell

    crystal = Crystal(
        cell_vectors=(primitive.cell * ureg("angstrom")),
        atom_r=(primitive.scaled_positions),
        atom_type=np.array(primitive.symbols),
        atom_mass=(primitive.masses * ureg("amu")),
    )

    # Phonopy s2p map uses supercell indices to identify atom;
    # e.g. p2s = [0, 4] s2p = [0, 0, 0, 0, 4, 4, 4, 4]
    # use argwhere to map back to primitive cell index
    # -> s2p n [0, 0, 0, 0, 1, 1, 1, 1]
    p2s_map = primitive.p2s_map
    s2p_map = np.argwhere(p2s_map[:, None] == primitive.s2p_map[None, :])[:, 0]

    fc_array, cell_origins = convert_fc_phases(
        phonon.force_constants,
        primitive.scaled_positions,
        supercell.scaled_positions @ phonon.supercell_matrix,
        p2s_map,
        s2p_map,
        phonon.supercell_matrix,
    )

    fc = ForceConstants(
        crystal=crystal,
        force_constants=(fc_array * ureg("eV angstrom^{-2}")),
        sc_matrix=phonon.supercell_matrix,
        cell_origins=cell_origins,
    )

    return fc


def main() -> None:
    args = get_parser().parse_args()
    phonon = get_phonopy(args.phonopy_yaml)
    fc = phonopy_to_euphonic(phonon)
    fc.to_json_file(args.output_json)


if __name__ == "__main__":
    main()
