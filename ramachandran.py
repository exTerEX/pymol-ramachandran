# Copyright 2021 Andreas Sagen Licensed under the Educational
# Community License, Version 2.0 (the "License"); you may not use this file
# except in compliance with the License. You may obtain a copy of the License at
#
# http://opensource.org/licenses/ECL-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

from pymol import cmd  # type: ignore
import matplotlib.pyplot as plt
import numpy

supported_amino_acids = (
    "arg", "his", "lys", "asp", "glu", "ser",
    "thr", "asn", "gln", "cys", "sec", "gly",
    "pro", "ala", "val", "ile", "leu", "met",
    "phe", "tyr", "trp"
)


def ramachandran(model: tuple = None, resn: tuple = None):
    """A simple function to generate a Ramachandran plot from a pymol instance. If model and/or resn is provided as an argument it limits the selection within the instance.

    Parameters
    ----------
    model : tuple, optional
        Model/enzyme within the pymol instance, by default None
    resn : tuple, optional
        Resn/residue within the model to limit phi-psi points, by default None

    Raises
    ------
    TypeError
        Trigger if model input isn't a string, list or tuple.
    TypeError
        Trigger if any model in a list or tuple isn't a string.
    ValueError
        Trigger if model isn't in the current pymol instance.
    ValueError
        Trigger if model isn't in the current pymol instance.
    TypeError
        Trigger if resn isn't a string, list or tuple.
    ValueError
        Trigger if resn isn't supported by pymol.
    TypeError
        Trigger if resn in a list or tuple isn't a string.
    ValueError
        Trigger if resn isn't supported by pymol.
    """

    # Validate that 'model' input is supported.
    if not isinstance(model, (str, list, tuple)) and model is not None:
        raise TypeError("'model' must be type 'str', 'list' or 'tuple'.")

    if isinstance(model, (list, tuple)):
        for element in model:
            if not isinstance(model, str):
                raise TypeError("'model' can only contain type 'str'.")

    local_model = cmd.get_names(type="public_objects")

    if isinstance(model, str) and not model in local_model:
        raise ValueError(f"Model '{model}' isn't in the pymol instance.")

    if isinstance(model, (list, tuple)) and not model in local_model:
        for element in model:
            if not element in local_model:
                raise ValueError(
                    f"Model '{element}' isn't in the pymol instance.")

        # Validate that 'resn' input is supported.
    if not isinstance(resn, (str, list, tuple)) and resn is not None:
        raise TypeError("'resn' must be type 'str', 'list' or 'tuple'.")

    if isinstance(resn, str) and not resn in supported_amino_acids:
        raise ValueError(f"'{resn}' isn't a supported amino acid residue.")

    if isinstance(resn, (list, tuple)):
        for residue in resn:
            if not isinstance(residue, str):
                raise TypeError("'resn' can only contain type 'str'.")

            if not residue in supported_amino_acids:
                raise ValueError(
                    f"'{residue}' isn't a supported amino acid residue.")

    if not model:
        model = local_model
        del local_model

    # Cycle through models
    for element in model:
        # Select model and/or from resn depending on condition.
        if resn is not None:
            cmd.select(f"sele_{element}", f"m. {element} & r. {resn}")
        else:
            cmd.select(f"sele_{element}", f"m. {element}")

        # Get all phi/psi values from entire structure
        phi_psi = cmd.phi_psi(selection=f"sele_{element}")

        # Remove intermediate selection
        cmd.delete(f"sele_{element}")

        # Seperate phi_psi into phi and psi.
        phi, psi = zip(*phi_psi.values())

        plt.figure(figsize=(8.5, 5.0), dpi=100)
        plt.title(f"Ramachandran plot (PDBid: {element.upper()})")

        plt.axes().set_aspect("equal", "box")

        plt.xlabel("\u03C6")
        plt.ylabel("\u03C8")

        min_bound = -180
        max_bound = +180

        plt.xlim(min_bound, max_bound)
        plt.ylim(min_bound, max_bound)

        ticks = numpy.arange(min_bound, max_bound + 1, 45, dtype=int)

        plt.xticks(ticks)
        plt.yticks(ticks)

        plt.axhline(y=0, color="k", lw=0.5)
        plt.axvline(x=0, color="k", lw=0.5)

        plt.grid(b=None, which="major",
                 axis="both", color="k", alpha=0.2)

        density_estimate_data = numpy.loadtxt(
            "data/density_estimate.csv", delimiter=",")

        density = numpy.log(numpy.rot90(density_estimate_data))
        plt.imshow(density, cmap="inferno", extent=(
            min_bound, max_bound) * 2, alpha=0.70)

        contour = numpy.rot90(numpy.fliplr(density_estimate_data))
        plt.contour(contour, colors="k", linewidths=0.5,
                    levels=[10**i for i in range(-7, 0)],
                    antialiased=True, extent=(min_bound, max_bound) * 2, alpha=0.55)

        if resn is not None:
            plt.scatter(phi, psi, marker=".", s=3, c="k")
        else:
            plt.scatter(phi, psi, marker=".", s=3, c="k")

        plt.show()


cmd.extend("ramachandran", ramachandran)
