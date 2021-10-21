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


"""Provide a sequencing sample model."""


from __future__ import annotations

import logging
from typing import Optional, Type

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from ..model import (
    CIGAROperation,
    Conflict,
    ConflictReliability,
    ConflictType,
    Plasmid,
    Sample,
    SampleReport,
)
from .sequence_alignment_service import SequenceAlignmentService


logger = logging.getLogger(__name__)


class SampleReportBuilder:
    def __init__(self, *, alignment_service: Type[SequenceAlignmentService], **kwargs):
        super().__init__(**kwargs)
        self.alignment_service = alignment_service
        self.report: Optional[SampleReport] = None

    def build(
        self,
        plasmid: Plasmid,
        sample: Sample,
        quality_threshold: float,
        window: int,
        **kwargs,
    ) -> SampleReport:
        self.report = SampleReport(
            sample=sample, quality_threshold=quality_threshold, window=window
        )
        if self._build_median():
            return self.report
        if self._build_trimmed():
            return self.report
        if self._build_alignment(plasmid, **kwargs):
            return self.report
        if self._build_conflicts():
            return self.report
        return self.report

    def _build_median(self) -> bool:
        self.report.median_quality = np.nanmedian(self.report.sample.phred_quality)
        if self.report.median_quality < self.report.quality_threshold:
            self.report.errors.append(
                f"The sample's ({self.report.sample.identifier}) median Phred quality "
                f"{self.report.median_quality} is below the threshold "
                f"{self.report.quality_threshold}."
            )
            return True
        return False

    def _build_trimmed(
        self,
        prefix: str = "",
        suffix: str = "_trimmed",
    ) -> bool:
        """Trim a sequence sample based on Phred quality."""
        index = np.arange(len(self.report.sample), dtype=int)
        self.report.smoothed = lowess(
            self.report.sample.phred_quality,
            index,
            frac=self.report.window / len(index),
            is_sorted=True,
            return_sorted=False,
        )
        mask = self.report.smoothed >= self.report.quality_threshold
        # Get min/max position that has a Phred quality higher than the threshold.

        self.report.trim_begin = index[mask][0]
        self.report.trim_end = index[mask][-1] + 1
        self.report.trimmed = Sample(
            identifier=f"{prefix}{self.report.sample.identifier}{suffix}",
            sequence=self.report.sample.sequence[
                self.report.trim_begin : self.report.trim_end
            ],
            phred_quality=self.report.sample.phred_quality[
                self.report.trim_begin : self.report.trim_end
            ],
        )
        return False

    def _build_alignment(self, plasmid: Plasmid, **kwargs) -> bool:
        self.report.alignment = self.alignment_service.align(
            query=self.report.trimmed.sequence, target=plasmid.sequence, **kwargs
        )
        self.report.strand = 1
        direction = "forward"
        note = (
            f"Alignment of trimmed sample {self.report.sample.identifier} to "
            f"plasmid {plasmid.identifier}."
        )
        if self.report.alignment.identity < 0.95:
            reverse = self.alignment_service.align(
                query=self.report.trimmed.sequence.reverse_complement(),
                target=plasmid.sequence,
                **kwargs,
            )
            if reverse.score > self.report.alignment.score:
                self.report.alignment = reverse
                self.report.strand = -1
                direction = "reverse"
                note = (
                    f"Alignment of reverse complement of trimmed sample "
                    f"{self.report.sample.identifier} to plasmid "
                    f"{plasmid.identifier}."
                )
        logger.info(
            "Add sample from %d to %d in %s direction.",
            self.report.alignment.target_begin,
            self.report.alignment.target_end,
            direction,
        )
        plasmid.add_feature(
            self.report.alignment.target_begin,
            self.report.alignment.target_end,
            type="misc_feature",
            strand=self.report.strand,
            qualifiers={
                "label": f"Alignment of {self.report.sample.identifier}",
                "note": note,
            },
        )
        return False

    def _build_conflicts(self) -> bool:
        target_index = self.report.alignment.target_begin
        query_index = self.report.alignment.query_begin
        events = list(self.report.alignment.iter_cigar())
        for idx, (span, code) in enumerate(events):
            if code is CIGAROperation.match:
                target_index += span
                query_index += span
            elif code is CIGAROperation.deletion:
                # Our logic requires that the preceding and following event are matches.
                assert events[idx - 1][1] is CIGAROperation.match
                assert events[idx + 1][1] is CIGAROperation.match
                quality = (
                    self.report.trimmed.phred_quality[query_index - 1]
                    + self.report.trimmed.phred_quality[query_index]
                ) / 2
                logger.info(
                    "Add deletion from %d to %d.", target_index, target_index + span
                )
                self.report.conflicts.append(
                    Conflict(
                        type=ConflictType.DELETION,
                        begin=target_index,
                        end=target_index + span,
                        reliability=ConflictReliability.HIGH
                        if quality >= self.report.quality_threshold
                        else ConflictReliability.LOW,
                        span=span,
                    )
                )
                target_index += span
            elif code is CIGAROperation.insertion:
                # Our logic requires that the preceding and following event are matches.
                assert events[idx - 1][1] is CIGAROperation.match
                assert events[idx + 1][1] is CIGAROperation.match
                quality = np.nanmean(
                    self.report.trimmed.phred_quality[query_index : query_index + span]
                )
                logger.info("Add insertion at %d.", target_index)
                self.report.conflicts.append(
                    Conflict(
                        type=ConflictType.INSERTION,
                        begin=target_index,
                        end=target_index,
                        reliability=ConflictReliability.HIGH
                        if quality >= self.report.quality_threshold
                        else ConflictReliability.LOW,
                        span=span,
                    )
                )
                query_index += span
            else:
                raise ValueError("Unhandled CIGAR event.")
        return False
