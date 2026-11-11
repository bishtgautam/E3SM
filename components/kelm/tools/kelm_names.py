"""Names that changed when ELMxx was renamed to kELM.

Output from a case made before the rename (compset I1850ELMXX / IELMXX,
longname ..._ELMXX_...) says 'elmxx' where a kELM case (I1850KELM / IKELM,
..._KELM_...) says 'kelm':

  files                 <case>.elmxx.h0.*, <case>.elmxx.r.*
  diagnostic labels     elmxx_in:, elmxx_out:, elmxxmap:, elmxxurb:, ...

The tools are written in kELM names. These helpers decide from a file's
COMPSET which name it uses and translate. The COMPSET comes from the case
that owns the file: the run dir's CASEROOT file, or an env_case.xml in the
file's directory or one of its parents. When no case can be found (a file
copied out of its run dir), both names are accepted -- a kELM file never
holds an 'elmxx' label, so mapping old to new is unambiguous.
"""

import os
import re

NEW = "kelm"
OLD = "elmxx"

_COMPSET_RE = re.compile(r'<entry id="COMPSET" value="([^"]*)"')


def _env_case(directory):
    """env_case.xml of the case that owns `directory`, or None."""
    d = os.path.abspath(directory)
    for _ in range(4):
        caseroot = os.path.join(d, "CASEROOT")  # written into every RUNDIR
        if os.path.isfile(caseroot):
            with open(caseroot) as f:
                candidate = os.path.join(f.read().strip(), "env_case.xml")
            if os.path.isfile(candidate):
                return candidate
        candidate = os.path.join(d, "env_case.xml")
        if os.path.isfile(candidate):
            return candidate
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return None


def compset(path):
    """COMPSET longname of the case that owns `path`, or None."""
    directory = path if os.path.isdir(path) else os.path.dirname(path) or "."
    env = _env_case(directory)
    if env is None:
        return None
    with open(env) as f:
        m = _COMPSET_RE.search(f.read())
    return m.group(1) if m else None


def component(path):
    """'elmxx' for a pre-rename case, 'kelm' for a kELM case, None if unknown."""
    cs = compset(path)
    if cs is None:
        return None
    return OLD if "_ELMXX" in cs else NEW


def label_mapper(path):
    """Function mapping a label as stored in `path` to its kELM name."""
    if component(path) == NEW:
        return lambda label: label

    def to_new(label):
        return NEW + label[len(OLD):] if label.startswith(OLD) else label
    return to_new


def file_names(path):
    """Component names to try in file globs for `path`, most likely first."""
    comp = component(path)
    return [comp] if comp else [NEW, OLD]
