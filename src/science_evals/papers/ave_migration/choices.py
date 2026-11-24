"""The fixed answer options of ave04 and the rule that picks the keyed one from the data.

Kept apart from reference.py so that the task specs and the tests can import them without
the analysis stack (scipy, scanpy).
"""
from __future__ import annotations

# Fixed option texts; the keyed option is chosen by decide_dkk1 from the measured AUCs.
AVE04_CHOICES = [
    "The data support DKK1 as a guidance cue: pooled over both stages, Dkk1 is clearly higher in the AVE than in the emVE, so the AVE releases DKK1 while it migrates.",
    "The data argue against AVE-derived DKK1 starting the migration: Dkk1 is no higher in the AVE than in the emVE at 5.5 dpc and becomes AVE-enriched only at 6.25 dpc.",
    "Dkk1 is enriched in the AVE at both stages about as strongly as Cer1, so it marks the AVE from its specification onwards and could guide the migration from the start.",
    "Dkk1 is AVE-enriched at 6.25 dpc and DKK1 is an established AVE guidance cue, so this enrichment confirms that it directs the migration that comes before this stage.",
    "Dkk1 is lower in the AVE than in the emVE at both stages, so any guidance by DKK1 would come from the emVE, pushing AVE cells away rather than being sent by them.",
    "Dkk1 is barely detected in the visceral endoderm at 5.5 dpc, so the data say nothing about its role at that stage, and the question has to wait for other experiments.",
]
NULL_BAND = 0.1  # an AUC within 0.5 +/- 0.1 counts as no enrichment
ENRICHED = 0.65  # an AUC at or above this counts as AVE-enriched


def decide_dkk1(aucs: dict, detected_e55: float) -> int:
    """The option (1-based) that the measured pattern supports."""
    d55, d625 = aucs["dkk1_auc_e55"], aucs["dkk1_auc_e625"]
    if detected_e55 < 0.2:
        return 6
    if abs(d55 - 0.5) <= NULL_BAND and d625 >= ENRICHED:
        return 2
    if d55 >= ENRICHED and d625 >= ENRICHED:
        return 3
    if d55 <= 1 - ENRICHED and d625 <= 1 - ENRICHED:
        return 5
    raise ValueError(f"no option matches the measured Dkk1 pattern {aucs}")
