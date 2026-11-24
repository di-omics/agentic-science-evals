"""One package per paper. Each exposes PAPER, DATA_FILES, build(data_dir) and,
optionally, anchors(data_dir). See docs/BUILDING.md for the contract."""
from __future__ import annotations

from dataclasses import dataclass, field

from ..controls import Control
from ..spec import Task

PAPERS = ("ave_migration", "cell_competition", "apoe_microglia")


@dataclass
class Built:
    """What a paper's build() returns for one task."""

    task: Task
    reference: dict  # the reference analysis's answer, computed from the data
    controls: list[Control]
    data_files: list[dict] = field(default_factory=list)  # DATA_FILES entries the agent receives
