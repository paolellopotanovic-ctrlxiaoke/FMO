import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from freeze_current_cycle_source import require_accepted_response


def test_only_structurally_accepted_label_free_response_can_replace_source():
    decision = {"action": "ACCEPT_LOCAL_REPAIR", "accepted": True,
                "affinity_labels_used": False, "scores_used": False, "checks": {"geometry": True}}
    require_accepted_response(decision)
    with pytest.raises(ValueError, match="label-free"):
        require_accepted_response({**decision, "scores_used": True})
    with pytest.raises(ValueError, match="label-free"):
        require_accepted_response({**decision, "checks": {"geometry": False}})
