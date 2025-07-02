import os
from re import (compile, match,)
import numpy as np
import pathlib
import pandas as pd
from typing import (Dict, List, Union)

H_TO_EV = 27.2114

pattern_force_eval = compile(r".*FORCE_EVAL\s+.*\s(-?\d+\.\d+)")

def search_energy_values(path: str) -> Dict[str, Dict[str, Union[float, Dict[str, float]]]]:
    """
    Scans the given directory and processes subfolders corresponding to different molecules. Each molecule folder 
    contains three subdirectories, each with a CP2K `.out` file representing the neutral, anion, and cation species.

    For each species, the script extracts energy values and constructs a nested dictionary of the following structure:
    
        {
            "MOLECULE_NAME": {
                "neutral": E_neutral,
                "anion": {
                    "vertical": E_anion_vertical,
                    "adiabatic": E_anion_adiabatic
                },
                "cation": {
                    "vertical": E_cation_vertical,
                    "adiabatic": E_cation_adiabatic
                }
            },
            ...
        }

    Parameters
    ----------
    path : str
        Path to the root directory containing molecule subdirectories.

    Returns
    -------
    Dict[str, Dict[str, Union[float, Dict[str, float]]]]
        A dictionary mapping molecule names to their corresponding energies for neutral, anionic, and cationic species.
    """

    if not isinstance(path, str):
        raise TypeError("Path must be a string.")
    elif str == "":
        raise ValueError("Path must be a non empty string.")

    main_folder_list = [f for f in os.scandir(path) if f.is_dir()]
    if not main_folder_list:
        raise IOError(f"No directories found in {path}.")

    energies = {}

    for folder in main_folder_list:
        subfolders_list = [s for s in os.scandir(folder) if s.is_dir()]
        subfolders_names_list = [s.name for s in subfolders_list]
        if not ("anion" in subfolders_names_list and "cation" in subfolders_names_list and "neutral" in subfolders_names_list):
            continue
        else:
            energies[folder.name] = {}
            for subfolder in subfolders_list:
                try:
                    out_file = [f.path for f in os.scandir(subfolder) if pathlib.Path(f).suffix == ".out" and pathlib.Path(f).stem != "job"][0]
                except:
                    FileNotFoundError(f"No valid .out file found in {subfolder.path}.")
                with open(out_file, "r") as infile:
                    energies_list = []
                    for line in infile:
                        found_energy = match(pattern_force_eval, line.strip())
                        if found_energy:
                            energies_list.append(float(found_energy.group(1)))
                    match subfolder.name:
                        case "neutral":
                            energies[folder.name][subfolder.name] = energies_list[-1]
                        case "anion"|"cation":
                            energies[folder.name][subfolder.name] = {}
                            energies[folder.name][subfolder.name]["vertical"] = energies_list[0]
                            energies[folder.name][subfolder.name]["adiabatic"] = energies_list[-1]
    return energies        

def relax_calculator(energies: Dict) -> None:
    """
    Computes relaxation energies and electronic properties from a nested energy dictionary 
    and exports the results to a CSV file.

    The function processes molecular energy data for neutral, cationic, and anionic species. 
    It calculates:
      - Vertical Ionization Potential (VIP)
      - Adiabatic Ionization Potential (AIP)
      - Vertical Electron Affinity (VEA)
      - Adiabatic Electron Affinity (AEA)
      - Relaxation energies for cations (lambda_plus)
      - Relaxation energies for anions (lambda_minus)

    All computed values are saved in a file named `report.csv`.

    Parameters
    ----------
    energies : Dict
        A nested dictionary of molecular energies with the following structure:
        {
            "Molecule1": {
                "neutral": float,
                "cation": {"vertical": float, "adiabatic": float},
                "anion": {"vertical": float, "adiabatic": float}
            },
            ...
        }

    Returns
    -------
    None
        The function saves the output as a CSV file and does not return any value.
    """

    if not isinstance(energies, dict):
        raise TypeError("The energies parameter must be a nested dictionary.")

    Molecules = []
    Neutral=[]
    Cation_vertical=[]
    Cation_adiabatic=[]
    Anion_vertical=[]
    Anion_adiabatic=[]

    for k,v in energies.items():
        Molecules.append(k)

        for k1,v1 in v.items():
            if k1 == "neutral":
                Neutral.append(v1)
            if k1 == "cation":
                for k2,v2 in v1.items():
                    if k2 == "vertical":
                        Cation_vertical.append(v2)
                    if k2 == "adiabatic":
                        Cation_adiabatic.append(v2)
            if k1 == "anion":
                for k2,v2 in v1.items():
                    if k2 == "vertical":
                        Anion_vertical.append(v2)
                    if k2 == "adiabatic":
                        Anion_adiabatic.append(v2)

    Neutral= np.array(Neutral)
    Cation_vertical = np.array(Cation_vertical)
    Cation_adiabatic = np.array(Cation_adiabatic)
    Anion_vertical = np.array(Anion_vertical)
    Anion_adiabatic = np.array(Anion_adiabatic) 
    VIP = (Cation_vertical - Neutral)*H_TO_EV
    AIP = (Cation_adiabatic - Neutral)*H_TO_EV
    VEA = (Neutral - Anion_vertical)*H_TO_EV
    AEA = (Neutral - Anion_adiabatic)*H_TO_EV
    lamb_plus = np.abs(VIP - AIP)
    lamb_minus = AEA - VEA
    
    df = pd.DataFrame({
        "Molecule": Molecules,
        "Neutral": Neutral,
        "Cation_vertical": Cation_vertical,
        "Cation_adiabatic": Cation_adiabatic,
        "Anion_vertical": Anion_vertical,
        "Anion_adiabatic": Anion_adiabatic,
        "VIP": VIP,
        "AIP": AIP,
        "VEA": VEA,
        "AEA": AEA,
        "lambda_plus": lamb_plus,
        "lambda_minus": lamb_minus
    })

    df.to_csv("report.csv", index=False, float_format="%.4f")