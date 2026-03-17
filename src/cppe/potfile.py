from __future__ import annotations

from pathlib import Path

from .pycppe import Multipole, Polarizability, Potential, multipole_components


_ANG2BOHR = 1.8897261246


class PotfileReader:
    def __init__(self, potfile_name: str):
        self._potfile = Path(potfile_name)
        if not self._potfile.is_file():
            raise RuntimeError("Potential file does not exist.")

    def read(self) -> list[Potential]:
        with self._potfile.open("r", encoding="utf-8") as handle:
            lines = handle.readlines()

        potentials: list[Potential] = []
        num_sites = 0
        i = 0

        while i < len(lines):
            line = lines[i]
            if "@COORDINATES" in line:
                i += 1
                num_sites = int(lines[i].split()[0])
                i += 1
                unit = lines[i].strip()
                if unit == "AA":
                    conversion = _ANG2BOHR
                elif unit == "AU":
                    conversion = 1.0
                else:
                    raise RuntimeError("Invalid unit for potential file.")

                for site_idx in range(num_sites):
                    i += 1
                    tokens = lines[i].split()
                    if len(tokens) < 4:
                        raise RuntimeError("Invalid coordinate line in potential file.")
                    element = tokens[0]
                    x = float(tokens[1]) * conversion
                    y = float(tokens[2]) * conversion
                    z = float(tokens[3]) * conversion
                    potentials.append(Potential(x, y, z, element, site_idx))

            elif "ORDER" in line:
                tokens = line.split()
                if len(tokens) == 2:
                    order = int(tokens[1])
                    i += 1
                    num_multipoles = int(lines[i].split()[0])
                    site_before = -1

                    for n_mul in range(num_multipoles):
                        i += 1
                        tokens = lines[i].split()
                        site_num = int(tokens[0]) - 1

                        if site_num != site_before + 1:
                            diff = site_num - site_before
                            for d in range(1, diff):
                                mul = Multipole(order)
                                for _ in range(multipole_components(order)):
                                    mul.add_value(0.0)
                                potentials[site_before + d].add_multipole(mul)

                        mul = Multipole(order)
                        expected = multipole_components(order)
                        if len(tokens) < expected + 1:
                            raise RuntimeError(
                                "Invalid multipole line in potential file."
                            )
                        for value_token in tokens[1 : expected + 1]:
                            mul.add_value(float(value_token))
                        mul.remove_trace()
                        potentials[site_num].add_multipole(mul)
                        site_before = site_num

                        if n_mul == num_multipoles - 1 and site_num != (num_sites - 1):
                            diff = num_sites - site_num
                            for d in range(1, diff):
                                missing_mul = Multipole(order)
                                for _ in range(multipole_components(order)):
                                    missing_mul.add_value(0.0)
                                potentials[site_num + d].add_multipole(missing_mul)

                elif len(tokens) == 3:
                    order_1 = int(tokens[1])
                    order_2 = int(tokens[2])
                    if order_1 != 1 or order_2 != 1:
                        raise RuntimeError(
                            "Only dipole-dipole polarizabilities are currently supported."
                        )
                    i += 1
                    num_polarizabilities = int(lines[i].split()[0])
                    for _ in range(num_polarizabilities):
                        i += 1
                        tokens = lines[i].split()
                        site_num = int(tokens[0]) - 1
                        expected = multipole_components(order_1 + order_2)
                        if len(tokens) < expected + 1:
                            raise RuntimeError(
                                "Invalid polarizability line in potential file."
                            )
                        values = [float(token) for token in tokens[1 : expected + 1]]
                        if potentials[site_num].is_polarizable:
                            raise RuntimeError(
                                "Potential can only have one polarizability."
                            )
                        potentials[site_num].set_polarizability(Polarizability(values))

                else:
                    raise RuntimeError("Invalid number in potfile ORDER.")

            elif "EXCLISTS" in line:
                i += 1
                num_exclusions = int(lines[i].split()[0])
                for _ in range(num_exclusions):
                    i += 1
                    tokens = lines[i].split()
                    if not tokens:
                        continue
                    site_num = int(tokens[0]) - 1
                    for token in tokens[1:]:
                        excl_site_number = int(token) - 1
                        if excl_site_number < 0:
                            continue
                        potentials[site_num].add_exclusion(excl_site_number)

            i += 1

        return potentials
