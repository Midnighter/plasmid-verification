# Copyright (c) 2021, Moritz E. Beber.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""Provide a sequencing sample report."""


from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Optional

import numpy as np

from .sample import Sample
from .sequence_alignment import SequenceAlignment


class ConflictType(str, Enum):
    """Define the possible types of alignment conflicts."""

    INSERTION = "insertion"
    DELETION = "deletion"
    CHANGE = "change"


class ConflictReliability(str, Enum):
    """Define the descriptors of the Phred quality in the sequence region."""

    HIGH = "high"
    LOW = "low"


class ConflictStatus(str, Enum):
    """Define the possible status types of an alignment conflict."""

    RESOLVED = "resolved"
    UNRESOLVED = "unresolved"
    POTENTIAL = "potential"


class EffectType(str, Enum):
    """Define the possible effects of sequence change."""

    FRAME_SHIFT = "frame shift"
    UKNOWN = "unknown"
    AA_CHANGE = "amino acid change"


@dataclass()
class Conflict:

    type: ConflictType
    begin: int
    end: int
    span: int
    reliability: ConflictReliability
    status: ConflictStatus = ConflictStatus.POTENTIAL
    num_confirmed: int = 0
    num_invalidated: int = 0

    def confirm(self) -> None:
        self.num_confirmed += 1

    def invalidate(self) -> None:
        self.num_invalidated += 1

    def evaluate_status(self) -> None:
        if self.num_confirmed > 0 and self.num_invalidated == 0:
            self.status = ConflictStatus.UNRESOLVED
        elif self.num_confirmed == 0 and self.num_invalidated > 0:
            self.status = ConflictStatus.RESOLVED
        else:
            self.status = ConflictStatus.POTENTIAL


class SampleReport:
    def __init__(
        self,
        *,
        sample: Sample,
        quality_threshold: float,
        window: int,
        trimmed: Optional[Sample] = None,
        smooth: Optional[np.ndarray] = None,
        alignment: Optional[SequenceAlignment] = None,
        strand: int = 1,
        median_quality: Optional[float] = None,
        conflicts: Iterable[Conflict] = (),
        errors: Iterable[str] = (),
        warnings: Iterable[str] = (),
        **kwargs
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self.sample = sample
        self.quality_threshold = quality_threshold
        self.window = window
        self.trimmed = trimmed
        self.smooth = smooth
        self.alignment: Optional[SequenceAlignment] = alignment
        self.strand = strand
        self.median_quality: Optional[float] = median_quality
        self.conflicts = list(conflicts)
        self.errors = list(errors)
        self.warnings = list(warnings)
