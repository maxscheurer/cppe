from __future__ import annotations

import cppe

from tests.cache import cache


def _clone_potential_with_shift(potential, index, shift):
    px = potential.x + shift[0]
    py = potential.y + shift[1]
    pz = potential.z + shift[2]
    cloned = cppe.Potential(px, py, pz, potential.element, index)

    for multipole in potential.multipoles:
        m = cppe.Multipole(multipole.k)
        for value in multipole.values:
            m.add_value(float(value))
        cloned.add_multipole(m)

    if potential.is_polarizable:
        cloned.set_polarizability(
            cppe.Polarizability(list(potential.polarizability.values))
        )

    return cloned


def replicated_waterbox_system(repeats=(4, 4, 4), spacing_bohr=18.0):
    base_potentials = cppe.cpp.PotfileReader("tests/potfiles/pna_6w.pot").read()
    base_molecule = cache.molecule["pna"]

    nx, ny, nz = repeats
    all_potentials = []
    all_atoms = cppe.Molecule()

    base_n_sites = len(base_potentials)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                shift = (ix * spacing_bohr, iy * spacing_bohr, iz * spacing_bohr)
                pot_offset = len(all_potentials)

                for i, pot in enumerate(base_potentials):
                    all_potentials.append(
                        _clone_potential_with_shift(pot, pot_offset + i, shift)
                    )

                for i, pot in enumerate(base_potentials):
                    for excl in pot.exclusions:
                        all_potentials[pot_offset + i].add_exclusion(pot_offset + excl)

                for atom in base_molecule:
                    all_atoms.append(
                        cppe.Atom(
                            atom.atomic_number,
                            atom.x + shift[0],
                            atom.y + shift[1],
                            atom.z + shift[2],
                        )
                    )

    n_polsites = len(cppe.get_polarizable_sites(all_potentials))
    label = f"waterbox_{nx}x{ny}x{nz}_{n_polsites}sites"
    return {
        "label": label,
        "molecule": all_atoms,
        "potentials": all_potentials,
        "n_polsites": n_polsites,
        "base_sites": base_n_sites,
    }
