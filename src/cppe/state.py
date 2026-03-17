from __future__ import annotations

import numpy as np

from . import cpp


_DEFAULT_OPTIONS = {
    "potfile": "potential.pot",
    "iso_pol": False,
    "induced_thresh": 1e-8,
    "maxiter": 50,
    "damp_induced": False,
    "damping_factor_induced": 2.1304,
    "damp_multipole": False,
    "damping_factor_multipole": 2.1304,
    "summation_induced_fields": "direct",
    "tree_ncrit": 64,
    "tree_expansion_order": 5,
    "theta": 0.5,
    "pe_border": False,
    "border_type": "remove",
    "border_rmin": 2.2,
    "border_nredist": -1,
    "border_redist_order": 1,
    "border_redist_pol": False,
}


def _polarizable_sites(potentials):
    return [p for p in potentials if p.is_polarizable]


class CppeState:
    def __init__(self, options=None, molecule=None, printer=None):
        if options is None:
            options = {}
        if molecule is None:
            molecule = cpp.Molecule()
        if printer is None:
            printer = lambda _: None

        unknown = set(options).difference(_DEFAULT_OPTIONS)
        if unknown:
            key = sorted(unknown)[0]
            raise ValueError(f"Option key '{key}' is invalid.")

        self._options = dict(_DEFAULT_OPTIONS)
        self._options.update(options)
        self._molecule = molecule
        self._printer = printer

        if self._options["pe_border"]:
            raise NotImplementedError(
                "Python backend currently does not support pe_border manipulation."
            )

        from . import PotfileReader

        self._potentials = list(PotfileReader(self._options["potfile"]).read())
        if self._options["iso_pol"]:
            for potential in self._potentials:
                if potential.is_polarizable:
                    potential.polarizability.make_isotropic()

        self._nuc_fields = np.zeros(3 * len(self._potentials), dtype=np.float64)
        self._multipole_fields = np.zeros(3 * len(self._potentials), dtype=np.float64)
        self._polsites = _polarizable_sites(self._potentials)
        self._induced_moments = np.zeros(3 * len(self._polsites), dtype=np.float64)
        self._make_guess = True

        self._energies = {
            "Electrostatic": {"Electronic": 0.0, "Nuclear": 0.0, "Multipoles": 0.0},
            "Polarization": {"Electronic": 0.0, "Nuclear": 0.0, "Multipoles": 0.0},
        }

    @property
    def options(self):
        return self._options

    @property
    def energies(self):
        return self._energies

    @property
    def potentials(self):
        return self._potentials

    @property
    def static_fields(self):
        return self._nuc_fields + self._multipole_fields

    @property
    def nuclear_fields(self):
        return self._nuc_fields

    @property
    def multipole_fields(self):
        return self._multipole_fields

    @property
    def total_energy(self):
        ele = self.energies["Electrostatic"]
        pol = self.energies["Polarization"]
        return sum(ele.values()) + sum(pol.values())

    @property
    def summary_string(self):
        return (
            "Polarizable Embedding Summary\n"
            f"Electrostatic total: {sum(self.energies['Electrostatic'].values()):.12f}\n"
            f"Polarization total: {sum(self.energies['Polarization'].values()):.12f}\n"
            f"Total energy: {self.total_energy:.12f}"
        )

    @property
    def positions(self):
        pos = np.zeros((len(self._potentials), 3), dtype=np.float64)
        for i, p in enumerate(self._potentials):
            pos[i, 0] = p.x
            pos[i, 1] = p.y
            pos[i, 2] = p.z
        return pos

    @property
    def positions_polarizable(self):
        pos = np.zeros((len(self._polsites), 3), dtype=np.float64)
        for i, p in enumerate(self._polsites):
            pos[i, 0] = p.x
            pos[i, 1] = p.y
            pos[i, 2] = p.z
        return pos

    def get_polarizable_site_number(self):
        return len(self._polsites)

    def set_potentials(self, potentials):
        self._potentials = list(potentials)
        self._polsites = _polarizable_sites(self._potentials)
        self._induced_moments = np.zeros(3 * len(self._polsites), dtype=np.float64)
        self._nuc_fields = np.zeros(3 * len(self._potentials), dtype=np.float64)
        self._multipole_fields = np.zeros(3 * len(self._potentials), dtype=np.float64)
        self._make_guess = True

    def calculate_static_energies_and_fields(self):
        from . import MultipoleExpansion, MultipoleFields, NuclearFields

        multipole_expansion = MultipoleExpansion(self._molecule, self._potentials)
        self.energies["Electrostatic"]["Nuclear"] = (
            multipole_expansion.interaction_energy()
        )
        self._nuc_fields = NuclearFields(self._molecule, self._potentials).compute()
        self._multipole_fields = MultipoleFields(
            self._potentials, self._options
        ).compute()

    def nuclear_interaction_energy_gradient(self):
        from . import MultipoleExpansion

        return MultipoleExpansion(self._molecule, self._potentials).nuclear_gradient()

    def nuclear_field_gradient(self):
        from . import NuclearFields

        return NuclearFields(self._molecule, self._potentials).nuclear_gradient()

    def update_induced_moments(self, elec_fields, elec_only=False):
        from . import InducedMoments

        elec_fields = np.asarray(elec_fields, dtype=np.float64)
        if elec_only:
            total_fields = elec_fields
        else:
            total_fields = elec_fields + self._nuc_fields + self._multipole_fields

        ind = InducedMoments(self._potentials, self._options)
        ind.set_print_callback(self._printer)
        self._induced_moments = ind.compute(
            total_fields, self._induced_moments, self._make_guess
        )

        if self._make_guess:
            self._make_guess = False

        epol_elec = -0.5 * float(np.dot(self._induced_moments, elec_fields))
        self.energies["Polarization"]["Electronic"] = epol_elec
        if not elec_only:
            epol_nuclear = -0.5 * float(np.dot(self._induced_moments, self._nuc_fields))
            epol_multipoles = -0.5 * float(
                np.dot(self._induced_moments, self._multipole_fields)
            )
            self.energies["Polarization"]["Nuclear"] = epol_nuclear
            self.energies["Polarization"]["Multipoles"] = epol_multipoles

    def induced_moments_eef(self):
        from . import InducedMoments

        n_pol = len(self._polsites)
        ret = np.zeros((3 * n_pol, 3), dtype=np.float64)
        fdn = np.zeros((3 * n_pol, 3), dtype=np.float64)
        for s in range(n_pol):
            l = 3 * s
            fdn[l + 0, 0] = 1.0
            fdn[l + 1, 1] = 1.0
            fdn[l + 2, 2] = 1.0

        ind = InducedMoments(self._potentials, self._options)
        ind.set_print_callback(self._printer)
        for a in range(3):
            ret[:, a] = ind.compute(fdn[:, a], True)
        return ret

    def get_induced_moments(self):
        return self._induced_moments.tolist()
