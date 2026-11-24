"""Fail-closed final-state audit for recovered early FMO monomer SCF failures."""
import math
import re


def final_fmo_convergence(deck, lines):
    groups = {}
    for name in ("CONTRL", "SCF", "FMOPRP", "FMO", "PCM"):
        found = re.findall(r"\$" + name + r"\s+(.*?)\$END", deck, re.S | re.I)
        if len(found) != 1:
            raise ValueError("one complete registered FMO2/DFTB input group required")
        groups[name] = found[0]
    if (not re.search(r"\bSCFTYP=RHF\b", groups["CONTRL"]) or
            not re.search(r"\bRUNTYP=ENERGY\b", groups["CONTRL"]) or
            not re.search(r"\bNBODY=2\b", groups["FMO"]) or
            not re.search(r"\bNGUESS=[23]\b", groups["FMOPRP"]) or
            not re.search(r"\bIFMO=-1\b", groups["PCM"]) or
            re.search(r"\bNLAYER=([2-9]|[1-9]\d+)", groups["FMO"]) or
            "$DFTB " not in deck):
        raise ValueError("early-failure audit supports registered single-layer RHF FMO2-DFTB/PCM only")
    fragment_match = re.search(r"\bNFRAG=(\d+)\b", groups["FMO"])
    scf_match = re.search(r"\bCONV=([0-9.EeDd+-]+)", groups["SCF"])
    fmo_match = re.search(r"\bCONV=([0-9.EeDd+-]+)", groups["FMOPRP"])
    if not all((fragment_match, scf_match, fmo_match)):
        raise ValueError("missing fragment count or convergence thresholds")
    count = int(fragment_match[1])
    scf_tolerance = float(scf_match[1].replace("D", "E").replace("d", "e"))
    fmo_tolerance = float(fmo_match[1].replace("D", "E").replace("d", "e"))
    if count < 2 or not 0 < scf_tolerance <= 1e-8 or not 0 < fmo_tolerance <= 1e-7:
        raise ValueError("unregistered or relaxed final thresholds")
    header = re.compile(r"RHF\s+(monomer|dimer)\s+(\d+)(?:\s+(\d+))?\s+corr=")
    iteration = re.compile(r"FMO iter\s*(\d+)\.(\d+),\s*(\d+) fragment\(s\) converged,.*maxE\s*([0-9.EeDd+-]+|\*+)\s*$")
    current = None
    cycle = {}
    final_cycle = None
    early_failures = []
    dimer_count = 0
    normal = False
    layer_converged = False
    final_summary = False
    free_monomer_guess = False
    overflow_iterations = []

    def finish(block):
        nonlocal dimer_count
        if block is None:
            return
        if block["kind"] == "monomer":
            if final_cycle is not None:
                raise ValueError("unexpected monomer after final SCC cycle")
            if block["fragment"] in cycle:
                raise ValueError("ambiguous repeated monomer within SCC cycle")
            cycle[block["fragment"]] = block
        else:
            if (final_cycle is None or block["failed"] or not block["converged"] or
                    block["energy"] is None or block["density_tolerance"] is None or
                    block["density_tolerance"] > 2 * scf_tolerance * 1.000001):
                raise ValueError("dimer has no complete strict final SCF convergence proof")
            dimer_count += 1

    for number, line in enumerate(lines, 1):
        text = line.strip()
        free_monomer_guess |= text == "Using free molecule guess for monomers."
        if ("TERMINATED -ABNORMALLY-" in text or "DIMER CALCULATION DID NOT CONVERGE" in text or
                "Aborting." in text or "FMO IS NOT CONVERGED" in text):
            raise ValueError("unrecoverable runtime, dimer or FMO convergence failure")
        if text.startswith("RHF"):
            match = header.match(text)
            if match:
                finish(current)
                current = {"kind": match[1], "fragment": int(match[2]), "failed": False,
                           "converged": False, "energy": None, "density_tolerance": None, "line": number}
                continue
        if "SCF IS UNCONVERGED" in text:
            if current is None or current["kind"] != "monomer" or final_cycle is not None:
                raise ValueError("failure outside recoverable initial monomer SCC")
            current["failed"] = True
            early_failures.append({"fragment": current["fragment"], "line": number})
        if current is not None:
            if text.startswith("DENSITY MATRIX CONV="):
                current["density_tolerance"] = float(text.split("=", 1)[1].replace("D", "E"))
            if text in {"DENSITY CONVERGED", "ENERGY CONVERGED", "DIIS CONVERGED"}:
                current["converged"] = True
            if text.startswith("FINAL RDFTB ENERGY IS"):
                value = float(text.split()[4])
                if (not math.isfinite(value) or value == 0) and not current["failed"]:
                    raise ValueError("invalid fragment SCF energy")
                current["energy"] = value if math.isfinite(value) and value != 0 else None
        if text.startswith("FMO iter"):
            match = iteration.match(text)
            if not match or int(match[2]) != 1 or final_cycle is not None:
                raise ValueError("unsupported or ambiguous SCC iteration")
            finish(current)
            current = None
            if set(cycle) != set(range(1, count + 1)):
                raise ValueError("incomplete monomer SCC cycle")
            overflow = match[4].startswith("*")
            if overflow:
                if int(match[3]) != 0 or final_cycle is not None:
                    raise ValueError("unresolved FMO energy residual cannot certify converged fragments")
                overflow_iterations.append(int(match[1]))
            if int(match[3]) == count:
                residual = float(match[4].replace("D", "E"))
                if not math.isfinite(residual) or abs(residual) > fmo_tolerance:
                    raise ValueError("final FMO energy residual exceeds requested tolerance")
                if any(block["failed"] or not block["converged"] or block["energy"] is None or
                       block["density_tolerance"] is None or
                       block["density_tolerance"] > 2 * scf_tolerance * 1.000001 for block in cycle.values()):
                    raise ValueError("not all final monomers strictly converged")
                final_cycle = int(match[1])
            cycle = {}
        if "Monomer SCF for layer" in text and "converged in" in text:
            layer_converged = "Monomer SCF for layer 1 converged in" in text
        if text.startswith("Layer 1 converged in"):
            final_summary = True
        normal |= "EXECUTION OF GAMESS TERMINATED NORMALLY" in text
    finish(current)
    if re.search(r"\bNGUESS=3\b", groups["FMOPRP"]) and not free_monomer_guess:
        raise ValueError("requested free-monomer initial guess is not verified in the log")
    if not (normal and layer_converged and final_summary and final_cycle is not None and dimer_count > 0):
        raise ValueError("missing full final monomer/dimer convergence and normal completion")
    return {"schema": "strict-final-fmo2-scf-convergence-v1", "n_fragments": count,
            "nonconverged_residual_overflow_iterations": overflow_iterations,
            "free_monomer_initial_guess_verified": free_monomer_guess,
            "final_scc_iteration": final_cycle, "n_strict_final_monomers": count,
            "n_strict_quantum_dimers": dimer_count, "early_monomer_failures": early_failures,
            "final_scf_requested_tolerance": scf_tolerance, "final_fmo_requested_tolerance": fmo_tolerance,
            "all_early_failures_replaced_by_strict_final_monomers": True}
