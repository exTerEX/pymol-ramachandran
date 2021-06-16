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
import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy
from typing import Union, NoReturn, Optional

supported_amino_acids = {
    "r": "arg", "h": "his", "k": "lys", "d": "asp", "e": "glu", "s": "ser",
    "t": "thr", "n": "asn", "q": "gln", "c": "cys", "u": "sec", "g": "gly",
    "p": "pro", "a": "ala", "v": "val", "i": "ile", "l": "leu", "m": "met",
    "f": "phe", "y": "tyr", "w": "trp"
}


def ramachandran(model: Optional[tuple] = None,
                 resn: Optional[Union[str, tuple]] = None,
                 cResn: Optional[bool] = False) -> NoReturn:
    """A simple function to generate a Ramachandran plot from a pymol instance. If model and/or resn is provided as an argument it limits the selection within the instance.

    Parameters
    ----------
    model : tuple, optional
        Model/enzyme within the pymol instance, by default None
    resn : tuple, optional
        Resn/residue within the model to limit phi-psi points, by default None
    cResn : bool, optional
        If each amino acid should have a unique color and label in legder.
    """

    local_model = cmd.get_names(type="public_objects")

    # -- Model input validation
    if model is not None:
        if isinstance(model, (str, list, tuple)):
            if isinstance(model, str):
                temp = []
                temp.append(model)
                model = temp
                del temp

            if isinstance(model, tuple):
                model = list(model)

            for value in model:
                if not isinstance(value, str):
                    raise TypeError("'model' can only contain type 'str'.")

                if value not in local_model:
                    cmd.fetch(value)
        else:
            raise TypeError("'model' must be type 'str', 'list' or 'tuple'.")

    # -- Residue input validation
    if resn is not None:
        # Residue input is either string, list or tuple, else throw a TypeError
        if isinstance(resn, (str, list, tuple)):
            # Convert string and tuple to list type
            if isinstance(resn, tuple):
                resn = list(resn)

            if isinstance(resn, str):
                temp = []
                temp.append(resn)
                resn = temp
                del temp

            # Ensure valid residue codes
            for key, value in enumerate(resn):
                value = value.lower()

                if not isinstance(value, str):
                    raise TypeError("'resn' can only contain type 'str'.")

                if len(value) == 3:
                    if value not in supported_amino_acids.values():
                        ValueError(
                            f"'{value}' isn't a supported amino acid residue.")
                elif len(value) == 1:
                    if value not in supported_amino_acids.keys():
                        ValueError(
                            f"'{value}' isn't a supported amino acid residue.")

                    resn[key] = supported_amino_acids[value]
                else:
                    ValueError(
                        f"'{value}' isn't a supported amino acid residue.")

        else:
            raise TypeError("'resn' must be type 'str', 'list' or 'tuple'.")

    # If model is None
    if not model:
        if not local_model:
            raise RuntimeError("Neither a choosen, nor a local model to use.")

        model = local_model
        del local_model

    # If residue is non, include all amino acid values
    if resn is None:
        resn = list(supported_amino_acids.values())

    sString = []
    sCount = []
    phi_psi = []

    # Cycle through models
    for element in model:
        # Select model and/or from resn depending on condition.
        for key, value in enumerate(resn):
            sString.append(f"sele_{element}_{value}")
            sCount.append(cmd.select(
                sString[key], f"m. {element} & r. {value}"))

        # Get all phi/psi values from entire structure
        for value in sString:
            phi_psi.append(cmd.phi_psi(selection=value))

        # Remove intermediate selection
        for value in sString:
            cmd.delete(value)

        # TODO: Move to lower-level fig + ax interaction with matplotlib

        # -- Plot
        fig, ax = matplotlib.pyplot.subplots(figsize=(8.5, 5.0), dpi=100)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.10, top=0.95)

        # Set title of plot
        ax.set_title(f"Ramachandran plot ({element.upper()})", fontsize=12)

        # Set same x and y scaling
        ax.set_aspect("equal", "box")

        # Set x- and y-label
        ax.set_xlabel("\u03C6", fontsize=12)
        # TODO: ylabel to horizontal
        ax.set_ylabel("\u03C8", fontsize=12)

        # Set x- and y-limit for plot
        min_bound = -180
        max_bound = +180

        ax.set_xlim(min_bound, max_bound)
        ax.set_ylim(min_bound, max_bound)

        # Specify x- and y-tick frequency
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(45))
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(n=5))

        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(45))
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(n=5))

        # Set ax-lines
        ax.axhline(y=0, color="k", lw=0.5)
        ax.axvline(x=0, color="k", lw=0.5)

        # Set grid
        ax.grid(color="k", alpha=0.2)

        # Load favorable density data
        favorable_density_data = numpy.loadtxt(
            "data/favorable_density.csv", delimiter=",")

        density = numpy.log(numpy.rot90(favorable_density_data))

        # Display favorable density as heat profile
        ax.imshow(density, cmap="inferno", extent=(
            min_bound, max_bound) * 2, alpha=0.70, interpolation="lanczos")

        # Plot contour lines
        contour = numpy.rot90(numpy.fliplr(favorable_density_data))
        plt.contour(
            contour, colors="k", linewidths=0.5,
            levels=[10 ** i for i in range(-7, 0)],
            antialiased=True, extent=(min_bound, max_bound) * 2, alpha=0.55)

        # Add phi-psi data
        for key, value in enumerate(phi_psi):
            try:
                phi, psi = zip(*value.values())
            except ValueError as E:
                if not sCount[key] == 0:
                    raise E
                else:
                    phi, psi = (None, None)

            # Apply independent color and ledger if specified
            if cResn:
                ax.scatter(
                    phi, psi, marker=".", s=3,
                    label=f"{resn[key][0].upper()}{resn[key][1:].lower()}")
                ax.legend(
                    bbox_to_anchor=(1.0, 0.9, 1.0, 0.125),
                    loc="upper left", borderaxespad=0.0, framealpha=0.0)
            else:
                plt.scatter(phi, psi, marker=".", s=3, c="k")

        fig.show()


cmd.extend("ramachandran", ramachandran)
