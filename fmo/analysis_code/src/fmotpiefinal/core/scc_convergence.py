"""Strict DFTB SCC controls and fragment-resolved failure evidence."""
import re


def scc_diis_deck(deck):
    groups = list(re.finditer(r"(?m)^\s*\$DFTB\s+.*?\$END", deck, re.S))
    if len(groups) != 1:
        raise ValueError("exactly one DFTB group required")
    group = groups[0]
    if re.search(r"\bITYPMX\s*=", group.group(), re.I):
        raise ValueError("explicit SCC mixer already registered")
    replacement = group.group().replace("$END", "ITYPMX=2 $END")
    if any(len(line) > 80 for line in replacement.splitlines()):
        raise ValueError("SCC control exceeds GAMESS line width")
    return deck[:group.start()] + replacement + deck[group.end():]


def mixer_evidence(text):
    markers = {"diis": "DIIS CHARGE UPDATE", "broyden": "BROYDEN'S CHARGE MIXING",
               "simple": "SIMPLE CHARGE MIXING"}
    found = [name for name, marker in markers.items() if marker in text]
    return {"mixers": found, "verified_scc_diis": found == ["diis"] and "UNKNOWN ITYPMX" not in text}


def fragment_failure_graph(deck, text):
    count = re.search(r"\bNFRAG\s*=\s*(\d+)", deck)
    block = re.search(r"\bFRGNAM\(1\)\s*=\s*(.*?)\bINDAT\(1\)", deck, re.S)
    if not count or not block:
        raise ValueError("named FMO fragment inventory required")
    names = [value.strip() for value in block.group(1).replace("\n", " ").split(",") if value.strip()]
    if len(names) != int(count.group(1)):
        raise ValueError("fragment inventory mismatch")
    current, events = None, []
    for line_number, line in enumerate(text.splitlines(), 1):
        header = re.search(r"^\s*(?:RHF|ROHF|UHF)\s+(monomer|dimer)\s+(\d+)(?:\s+(\d+))?", line)
        if header:
            indices = [int(header.group(2))]
            if header.group(1) == "dimer":
                if header.group(3) is None:
                    raise ValueError("missing dimer partner")
                indices.append(int(header.group(3)))
            if any(index < 1 or index > len(names) for index in indices):
                raise ValueError("fragment header outside inventory")
            current = {"kind": header.group(1), "indices": indices,
                       "names": [names[index - 1] for index in indices], "header_line": line_number}
        if "SCF IS UNCONVERGED" in line:
            events.append({"line": line_number, "last_explicit_header": current})
    return {"fragment_nodes": names, "unconverged_events": events,
            "header_context_not_unique_causal_attribution": True, "mixer": mixer_evidence(text)}
