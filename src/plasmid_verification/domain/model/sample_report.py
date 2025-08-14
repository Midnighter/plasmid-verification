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


class StrandDirection(int, Enum):
    """Define the possible strand directions of a plasmid feature."""

    FORWARD = 1
    REVERSE = -1
    UNKNOWN = 0


class ConflictType(str, Enum):
    """Define the possible types of alignment conflicts."""

    INSERTION = "insertion"
    DELETION = "deletion"
    VARIANT = "variant"


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
    plasmid_begin: int
    plasmid_end: int
    sample_begin: int
    sample_end: int
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
        # TODO: call unresolved point mutation or SNP
        if self.num_confirmed > 0 and self.num_invalidated == 0:
            # TODO: Also include single read conflicts of high quality
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
        trimmed: Optional[Sample] = None,
        trim_begin: Optional[int] = None,
        trim_end: Optional[int] = None,
        smoothed: Optional[np.ndarray] = None,
        alignment: Optional[SequenceAlignment] = None,
        strand: StrandDirection = StrandDirection.FORWARD,
        median_quality: Optional[float] = None,
        conflicts: Iterable[Conflict] = (),
        errors: Iterable[str] = (),
        warnings: Iterable[str] = (),
        trim_kwargs: Optional[dict] = None,
        alignment_kwargs: Optional[dict] = None,
        **kwargs,
    ) -> None:
        """"""
        super().__init__(**kwargs)
        self.sample = sample
        self.quality_threshold = quality_threshold
        self.trimmed = trimmed
        self.trim_begin = trim_begin
        self.trim_end = trim_end
        self.smoothed = smoothed
        self.alignment: Optional[SequenceAlignment] = alignment
        self.strand = strand
        self.median_quality: Optional[float] = median_quality
        self.conflicts = list(conflicts)
        self.errors = list(errors)
        self.warnings = list(warnings)
        self.trim_kwargs = {} if trim_kwargs is None else trim_kwargs
        self.alignment_kwargs = {} if alignment_kwargs is None else alignment_kwargs

    def judge_reliability(self, score: float) -> ConflictReliability:
        return (
            ConflictReliability.HIGH
            if score >= self.quality_threshold
            else ConflictReliability.LOW
        )
